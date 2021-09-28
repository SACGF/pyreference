#!/usr/bin/env python

from __future__ import print_function, absolute_import

import HTSeq
import abc
import gzip
import json
import os
from argparse import ArgumentParser
from collections import defaultdict, Counter

from pyreference.settings import CHROM, START, END, STRAND, IS_CODING, \
    PYREFERENCE_JSON_VERSION_KEY, PYREFERENCE_JSON_VERSION
from pyreference.utils.file_utils import name_from_file_name, file_md5sum


class GFFParser(abc.ABC):
    FEATURES = {}  # Set of features that we process, override in subclasses

    def __init__(self, filename, discard_contigs_with_underscores=True):
        self.filename = filename
        self.discard_contigs_with_underscores = discard_contigs_with_underscores

        self.genes_by_id = {}
        self.transcripts_by_id = {}
        self.gene_id_by_name = {}
        self.gene_ids_by_biotype = defaultdict(set)
        # Store CDS in separate dict as we don't need to write as JSON
        self.transcript_cds_by_id = {}

    @abc.abstractmethod
    def handle_feature(self, feature):
        pass

    def parse(self):
        for feature in HTSeq.GFF_Reader(self.filename):
            if feature.type not in self.FEATURES:
                continue

            try:
                chrom = feature.iv.chrom
                if self.discard_contigs_with_underscores and not chrom.startswith("NC_") and "_" in chrom:
                    discarded_contigs[chrom] += 1
                    continue
                self.handle_feature(feature)
            except Exception as e:
                print("Could not parse '%s': %s" % (feature, e))
                raise e

    def finish(self):
        self._add_coding_and_utr_features()

    @staticmethod
    def _create_gene(gene_name, feature):
        return {
            "name": gene_name,
            "transcripts": set(),
            "biotype": set(),
            CHROM: feature.iv.chrom,
            START: feature.iv.start,
            END: feature.iv.end,
            STRAND: feature.iv.strand
        }

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
    def _get_biotype_from_transcript_id(transcript_id):
        biotypes_by_transcript_id_start = {"NM_": "protein_coding", "NR_": "non_coding"}
        for (start, biotype) in biotypes_by_transcript_id_start.items():
            if transcript_id.startswith(start):
                return biotype

        if "tRNA" in transcript_id:
            return "tRNA"
        return "N/A"

    @staticmethod
    def _update_extents(genomic_region_dict, feature):
        start = genomic_region_dict[START]
        if feature.iv.start < start:
            genomic_region_dict[START] = feature.iv.start

        end = genomic_region_dict[END]
        if feature.iv.end > end:
            genomic_region_dict[END] = feature.iv.end

    def _add_coding_and_utr_features(self):
        """ Add 5PUTR/3PUTR features to coding transcripts """

        for transcript_id, transcript in self.transcripts_by_id.items():
            cds_extent = self.transcript_cds_by_id.get(transcript_id)
            if cds_extent:
                transcript[IS_CODING] = 1
                features_by_type = transcript["features_by_type"]

                (left, right) = ("5PUTR", "3PUTR")
                if transcript[STRAND] == '-':  # Switch
                    (left, right) = (right, left)

                cds_start = cds_extent[START]
                cds_end = cds_extent[END]

                for exon in features_by_type["exon"]:
                    exon_start = exon[START]
                    exon_end = exon[END]

                    if exon_start < cds_start:
                        end_non_coding = min(cds_start, exon_end)
                        utr_feature = {START: exon_start,
                                       END: end_non_coding}
                        features_by_type[left].append(utr_feature)

                    if exon_end > cds_end:
                        start_non_coding = max(cds_end, exon_start)
                        utr_feature = {START: start_non_coding,
                                       END: exon_end}
                        features_by_type[right].append(utr_feature)

    def get_data(self):
        self.parse()
        self.finish()

        return {
            PYREFERENCE_JSON_VERSION_KEY: PYREFERENCE_JSON_VERSION,
            "reference_gtf": {"path": os.path.abspath(self.filename),
                              "md5sum": file_md5sum(self.filename)},
            "genes_by_id": self.genes_by_id,
            "transcripts_by_id": self.transcripts_by_id,
            "gene_id_by_name": self.gene_id_by_name,
            "gene_ids_by_biotype": self.gene_ids_by_biotype,
        }


class GTFParser(GFFParser):
    """ GTF (GFF2) - used by Ensembl, @see http://gmod.org/wiki/GFF2

        GFF2 only has 2 levels of feature hierarchy, so we have to build or 3 levels of gene/transcript/exons ourselves
    """
    CODING_FEATURES = {"CDS", "start_codon", "stop_codon"}
    FEATURES = CODING_FEATURES | {"exon"}

    def __init__(self, *args, **kwargs):
        super(GTFParser, self).__init__(*args, **kwargs)

    def handle_feature(self, feature):
        gene_id = feature.attr["gene_id"]
        gene_name = feature.attr.get("gene_name")  # Non mandatory - Ensembl doesn't have on some RNAs
        if gene_name:
            self.gene_id_by_name[gene_name] = gene_id  # TODO: Check for dupes?

        gene = self.genes_by_id.get(gene_id)
        if gene is None:
            gene = self._create_gene(gene_name, feature)
            self.genes_by_id[gene_id] = gene
        else:
            self._update_extents(gene, feature)

        transcript_id = feature.attr["transcript_id"]
        gene["transcripts"].add(transcript_id)
        transcript = self.transcripts_by_id.get(transcript_id)
        if transcript is None:
            transcript = self._create_transcript(feature)
            self.transcripts_by_id[transcript_id] = transcript
        else:
            self._update_extents(transcript, feature)

        # No need to store chrom/strand for each feature, will use transcript
        feature_dict = {START: feature.iv.start,
                        END: feature.iv.end, }

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

        biotype = feature.attr.get("gene_biotype")
        if biotype is None:
            biotype = self._get_biotype_from_transcript_id(transcript_id)

        if biotype:
            gene["biotype"].add(biotype)
            transcript["biotype"].add(biotype)

        self.gene_ids_by_biotype[biotype].add(gene_id)


class GFF3Parser(GFFParser):
    """ GFF3 - Used by RefSeq, @see https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md"""

    GFF3_GENES = {"gene", "pseudogene"}
    GFF3_TRANSCRIPTS = {"mRNA", "ncRNA"}
    GFF3_TRANSCRIPTS_DATA = {"exon", "CDS", "cDNA_match"}
    FEATURES = GFF3_GENES | GFF3_TRANSCRIPTS | GFF3_TRANSCRIPTS_DATA

    def __init__(self, *args, **kwargs):
        super(GFF3Parser, self).__init__(*args, **kwargs)

    def handle_feature(self, feature):
        pass


def handle_args():
    parser = ArgumentParser(description='Build a json.gz file for pyreference')
    parser.add_argument("--discard-contigs-with-underscores", action='store_true', default=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--gtf', help='GTF (Gene Transfer Format) filename')
    group.add_argument('--gff3', help='GFF3 (Gene Feature Format) filename')
    return parser.parse_args()


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
    discarded_contigs = Counter()

    parser = parser_factory(args.gtf, args.gff3,
                            discard_contigs_with_underscores=args.discard_contigs_with_underscores)
    data = parser.get_data()

    genes_json_gz = name_from_file_name(parser.filename) + ".json.gz"
    with gzip.open(genes_json_gz, 'w') as outfile:
        json_str = json.dumps(data, cls=SortedSetEncoder, sort_keys=True)  # Sort so diffs work
        outfile.write(json_str.encode('ascii'))

    print("Wrote:", genes_json_gz)
    if discarded_contigs:
        print("Discarded contigs: %s" % discarded_contigs)


if __name__ == '__main__':
    main()
