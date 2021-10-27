#!/usr/bin/env python

from __future__ import print_function, absolute_import

import HTSeq
import abc
import gzip
import json
import logging
import os
from argparse import ArgumentParser
from collections import defaultdict, Counter

from pyreference.settings import CHROM, START, END, STRAND, IS_CODING, \
    PYREFERENCE_JSON_VERSION_KEY, PYREFERENCE_JSON_VERSION
from pyreference.utils.file_utils import name_from_file_name, file_md5sum


class GFFParser(abc.ABC):
    CODING_FEATURES = {"CDS", "start_codon", "stop_codon"}
    FEATURE_ALLOW_LIST = {}
    FEATURE_IGNORE_LIST = {"biological_region", "chromosome", "region", "scaffold", "supercontig"}

    def __init__(self, filename, discard_contigs_with_underscores=True):
        self.filename = filename
        self.discard_contigs_with_underscores = discard_contigs_with_underscores

        self.discarded_contigs = Counter()
        self.genes_by_id = {}
        self.transcripts_by_id = {}
        self.gene_id_by_name = {}
        # Store CDS in separate dict as we don't need to write as JSON
        self.transcript_cds_by_id = {}

    @abc.abstractmethod
    def handle_feature(self, feature):
        pass

    def parse(self):
        for feature in HTSeq.GFF_Reader(self.filename):
            if self.FEATURE_ALLOW_LIST and feature.type not in self.FEATURE_ALLOW_LIST:
                continue
            if feature.type in self.FEATURE_IGNORE_LIST:
                continue

            try:
                chrom = feature.iv.chrom
                if self.discard_contigs_with_underscores and not chrom.startswith("NC_") and "_" in chrom:
                    self.discarded_contigs[chrom] += 1
                    continue
                self.handle_feature(feature)
            except Exception as e:
                print("Could not parse '%s': %s" % (feature.get_gff_line(), e))
                raise e

    def finish(self):
        self._add_coding_and_utr_features()

        if self.discarded_contigs:
            print("Discarded contigs: %s" % self.discarded_contigs)


    @staticmethod
    def _create_gene(gene_name, feature):
        biotypes = set()

        gene = {
            "name": gene_name,
            "transcripts": set(),
            "biotype": biotypes,
            CHROM: feature.iv.chrom,
            START: feature.iv.start,
            END: feature.iv.end,
            STRAND: feature.iv.strand
        }

        # Attempt to get some biotypes in there if available
        if feature.type == "gene":
            gene_version = feature.attr.get("version")
            biotype = feature.attr.get("biotype")
            description = feature.attr.get("description")
            if description:
                gene["description"] = description
        else:
            gene_version = feature.attr.get("gene_version")
            biotype = feature.attr.get("gene_biotype")

        if biotype:
            biotypes.add(biotype)

        if gene_version:
            gene["version"] = int(gene_version)
        return gene

    @staticmethod
    def _create_transcript(feature):
        return {
            "features_by_type": defaultdict(list),
            "biotype": set(),
            CHROM: feature.iv.chrom,
            START: feature.iv.start,
            END: feature.iv.end,
            STRAND: feature.iv.strand,
            IS_CODING: 0
        }

    @staticmethod
    def _store_other_chrom(data, feature):
        other_chroms = data.get("other_chroms", set())
        other_chroms.add(feature.iv.chrom)
        data["other_chroms"] = other_chroms

    @staticmethod
    def _get_biotype_from_transcript_id(transcript_id):
        biotypes_by_transcript_id_start = {"NM_": "protein_coding", "NR_": "non_coding"}
        for (start, biotype) in biotypes_by_transcript_id_start.items():
            if transcript_id.startswith(start):
                return biotype

        if "tRNA" in transcript_id:
            return "tRNA"
        return "N/A"

    def _add_transcript_data(self, transcript_id, transcript, feature):
        if feature.iv.chrom != transcript[CHROM]:
            self._store_other_chrom(transcript, feature)
            return

        feature_dict = {START: feature.iv.start,
                        END: feature.iv.end}
        if feature.type == "cDNA_match":
            target = feature.attr.get("Target")
            t_cols = target.split()
            feature_dict["cdna_start"] = int(t_cols[1])
            feature_dict["cdna_end"] = int(t_cols[2])
            if len(t_cols) == 4 and t_cols[3] != '+':  # Default is '+', so only store '-'
                feature_dict["cdna_strand"] = t_cols[3]
            gap = feature.attr.get("Gap")
            if gap:
                feature_dict["gap"] = gap

        transcript["features_by_type"][feature.type].append(feature_dict)
        if feature.type in self.CODING_FEATURES:
            cds_extent = self.transcript_cds_by_id.get(transcript_id)
            if cds_extent is None:
                cds_extent = {START: feature.iv.start,
                              END: feature.iv.end}
                self.transcript_cds_by_id[transcript_id] = cds_extent
            else:
                cds_extent[START] = min(cds_extent[START], feature.iv.start)
                cds_extent[END] = max(cds_extent[END], feature.iv.end)

    def _add_coding_and_utr_features(self):
        """ Add 5PUTR/3PUTR features to coding transcripts

            Ensembl GTFs have 'five_prime_UTR' features (similar to CDS etc) but we make this for GFFs that
            don't have those features
        """

        for transcript_id, transcript in self.transcripts_by_id.items():
            cds_extent = self.transcript_cds_by_id.get(transcript_id)
            if cds_extent:
                transcript[IS_CODING] = 1
                features_by_type = transcript["features_by_type"]

                (left, right) = ("5PUTR", "3PUTR")
                if transcript[STRAND] == '-':  # Switch
                    (left, right) = (right, left)

                cds_min = cds_extent[START]
                cds_max = cds_extent[END]

                transcript["cds_start"] = cds_min
                transcript["cds_end"] = cds_max

                # exon is in stranded order
                for exon in features_by_type["exon"]:
                    exon_start = exon[START]
                    exon_end = exon[END]

                    if exon_start < cds_min:
                        end_non_coding = min(cds_min, exon_end)
                        utr_feature = {START: exon_start,
                                       END: end_non_coding}
                        features_by_type[left].append(utr_feature)

                    if exon_end > cds_max:
                        start_non_coding = max(cds_max, exon_start)
                        utr_feature = {START: start_non_coding,
                                       END: exon_end}
                        features_by_type[right].append(utr_feature)

    def get_data(self):
        self.parse()
        self.finish()

        gene_ids_by_biotype = defaultdict(set)
        for gene_id, gene in self.genes_by_id.items():
            for biotype in gene["biotype"]:
                gene_ids_by_biotype[biotype].add(gene_id)

        return {
            PYREFERENCE_JSON_VERSION_KEY: PYREFERENCE_JSON_VERSION,
            "reference_gtf": {"path": os.path.abspath(self.filename),
                              "md5sum": file_md5sum(self.filename)},
            "genes_by_id": self.genes_by_id,
            "transcripts_by_id": self.transcripts_by_id,
            "gene_id_by_name": self.gene_id_by_name,
            "gene_ids_by_biotype": gene_ids_by_biotype,
        }


class GTFParser(GFFParser):
    """ GTF (GFF2) - used by Ensembl, @see http://gmod.org/wiki/GFF2

        GFF2 only has 2 levels of feature hierarchy, so we have to build or 3 levels of gene/transcript/exons ourselves
    """
    GTF_TRANSCRIPTS_DATA = GFFParser.CODING_FEATURES | {"exon"}
    FEATURE_ALLOW_LIST = GTF_TRANSCRIPTS_DATA | {"gene"}

    def __init__(self, *args, **kwargs):
        super(GTFParser, self).__init__(*args, **kwargs)

    def handle_feature(self, feature):
        gene_id = feature.attr["gene_id"]
        # Non mandatory - Ensembl doesn't have on some RNAs
        gene_name = None
        if feature.type == "gene":
            gene_name = feature.attr.get("Name")
        else:
            gene_name = feature.attr.get("gene_name")
        if gene_name:
            self.gene_id_by_name[gene_name] = gene_id  # Shouldn't be dupes per file

        gene = self.genes_by_id.get(gene_id)
        if gene is None:
            gene = self._create_gene(gene_name, feature)
            self.genes_by_id[gene_id] = gene
        else:
            self._update_extents(gene, feature)

        transcript_id = feature.attr.get("transcript_id")
        transcript_version = feature.attr.get("transcript_version")
        if transcript_version:
            transcript_id += "." + transcript_version

        if transcript_id:
            gene["transcripts"].add(transcript_id)
            transcript = self.transcripts_by_id.get(transcript_id)
            if transcript is None:
                transcript = self._create_transcript(feature)
                self.transcripts_by_id[transcript_id] = transcript
            else:
                self._update_extents(transcript, feature)

            # No need to store chrom/strand for each feature, will use transcript
            if feature.type in self.GTF_TRANSCRIPTS_DATA:
                self._add_transcript_data(transcript_id, transcript, feature)

            biotype = feature.attr.get("gene_biotype")
            if biotype is None:
                biotype = self._get_biotype_from_transcript_id(transcript_id)

            if biotype:
                gene["biotype"].add(biotype)
                transcript["biotype"].add(biotype)

    @staticmethod
    def _update_extents(genomic_region_dict, feature):
        if feature.iv.chrom == genomic_region_dict[CHROM]:
            start = genomic_region_dict[START]
            if feature.iv.start < start:
                genomic_region_dict[START] = feature.iv.start

            end = genomic_region_dict[END]
            if feature.iv.end > end:
                genomic_region_dict[END] = feature.iv.end
        else:
            self._store_other_chrom(genomic_region_dict, feature)


class GFF3Parser(GFFParser):
    """ GFF3 - Used by RefSeq, @see https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

        GFF3 support arbitrary hierarchy

    """

    GFF3_GENES = {"gene", "pseudogene"}
    GFF3_TRANSCRIPTS_DATA = {"exon", "CDS", "cDNA_match", "five_prime_UTR", "three_prime_UTR"}

    def __init__(self, *args, **kwargs):
        super(GFF3Parser, self).__init__(*args, **kwargs)
        self.gene_id_by_feature_id = defaultdict()
        self.transcript_id_by_feature_id = defaultdict()

    def handle_feature(self, feature):
        parent_id = feature.attr.get("Parent")
        # Genes never have parents
        # RefSeq genes are always one of GFF3_GENES, Ensembl has lots of different types (lincRNA_gene etc)
        # Ensembl treats pseudogene as a transcript (has parent)
        if parent_id is None and (feature.type in self.GFF3_GENES or "gene_id" in feature.attr):
            gene_id = feature.attr.get("gene_id")
            dbxref = self._get_dbxref(feature)
            if not gene_id:
                gene_id = dbxref.get("GeneID")
            if not gene_id:
                raise ValueError("Could not obtain 'gene_id', tried 'gene_id' and 'Dbxref[GeneID]'")

            gene_name = feature.attr.get("Name")
            # Gene can have multiple loci, thus entries in GFF, keep original so all transcripts are added
            gene = self.genes_by_id.get(gene_id)
            if gene is None:
                gene = self._create_gene(gene_name, feature)
                # If a gene already exists - then need to merge it...
                self.genes_by_id[gene_id] = gene

            hgnc = dbxref.get("HGNC")
            if hgnc:
                gene["HGNC"] = hgnc

            if gene_name:
                self.gene_id_by_name[gene_name] = gene_id
            self.gene_id_by_feature_id[feature.attr["ID"]] = gene_id
        else:
            if feature.type in self.GFF3_TRANSCRIPTS_DATA:
                if feature.type == 'cDNA_match':
                    target = feature.attr["Target"]
                    transcript_id = target.split()[0]
                else:
                    # Some exons etc may be for miRNAs that have no transcript ID, so skip those (won't have parent)
                    if parent_id:
                        transcript_id = self.transcript_id_by_feature_id.get(parent_id)
                    else:
                        logging.warning("Transcript data has no parent: %s" % feature.get_gff_line())
                        transcript_id = None

                if transcript_id:
                    transcript = self.transcripts_by_id[transcript_id]
                    self._handle_transcript_data(transcript_id, transcript, feature)
            else:
                # There are so many different transcript ontology terms just taking everything that
                # has a transcript_id and is child of gene (ie skip miRNA etc that is child of primary_transcript)
                transcript_id = feature.attr.get("transcript_id")
                if transcript_id:
                    transcript_version = feature.attr.get("version")
                    if transcript_version:
                        transcript_id += "." + transcript_version
                    assert parent_id is not None
                    gene_id = self.gene_id_by_feature_id.get(parent_id)
                    if not gene_id:
                        raise ValueError("Don't know how to handle feature type %s (not child of gene)" % feature.type)
                    gene = self.genes_by_id[gene_id]
                    self._handle_transcript(gene, transcript_id, feature)

    @staticmethod
    def _get_dbxref(feature):
        """ RefSeq stores attribute with more keys, eg: 'Dbxref=GeneID:7840,HGNC:HGNC:428,MIM:606844' """
        dbxref = {}
        dbxref_str = feature.attr.get("Dbxref")
        if dbxref_str:
            dbxref = dict(d.split(":", 1) for d in dbxref_str.split(","))
        return dbxref

    def _handle_transcript(self, gene, transcript_id, feature):
        """ Sometimes we can get multiple transcripts in the same file - just taking 1st """
        if transcript_id not in self.transcripts_by_id:
            # print("_handle_transcript(%s, %s)" % (gene, feature))
            transcript = self._create_transcript(feature)
            biotype = self._get_biotype_from_transcript_id(transcript_id)
            gene["transcripts"].add(transcript_id)
            gene["biotype"].add(biotype)
            transcript["biotype"].add(biotype)
            partial = feature.attr.get("partial")
            if partial:
                transcript["partial"] = 1
            self.transcripts_by_id[transcript_id] = transcript
        self.transcript_id_by_feature_id[feature.attr["ID"]] = transcript_id

    def _handle_transcript_data(self, transcript_id, transcript, feature):
        self._add_transcript_data(transcript_id, transcript, feature)


def handle_args():
    parser = ArgumentParser(description='Build a json.gz file for pyreference')
    parser.add_argument("--discard-contigs-with-underscores", action='store_true', default=True)
    parser.add_argument('--url', help='URL (source of GFF) to store in "reference_gtf.url"')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--gtf', help='GTF (Gene Transfer Format) filename')
    group.add_argument('--gff3', help='GFF3 (Gene Feature Format) filename')
    args = parser.parse_args()
    if not (args.gtf or args.gff3):
        parser.error("You must specify either --gtf or --gff3")
    return args


def parser_factory(gtf=None, gff3=None, discard_contigs_with_underscores=True):
    if gtf:
        parser = GTFParser(gtf, discard_contigs_with_underscores)
    else:
        parser = GFF3Parser(gff3, discard_contigs_with_underscores)
    return parser


class SortedSetEncoder(json.JSONEncoder):
    """ Dump set as list, from: https://stackoverflow.com/a/8230505/295724 """

    def default(self, obj):
        if isinstance(obj, set):
            return list(sorted(obj))
        return json.JSONEncoder.default(self, obj)


def main():
    args = handle_args()
    parser = parser_factory(args.gtf, args.gff3,
                            discard_contigs_with_underscores=args.discard_contigs_with_underscores)
    data = parser.get_data()
    if args.url:
        data["reference_gtf"]["url"] = args.url

    genes_json_gz = name_from_file_name(parser.filename) + ".json.gz"
    with gzip.open(genes_json_gz, 'w') as outfile:
        json_str = json.dumps(data, cls=SortedSetEncoder, sort_keys=True)  # Sort so diffs work
        outfile.write(json_str.encode('ascii'))

    print("Wrote:", genes_json_gz)


if __name__ == '__main__':
    main()
