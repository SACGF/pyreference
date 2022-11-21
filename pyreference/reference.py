from __future__ import print_function, absolute_import

import HTSeq
from bioutils.assemblies import make_ac_name_map
from collections import defaultdict
from deprecation import deprecated
from functools import reduce
import gzip
import json
from lazy import lazy
import operator
import os
from pyreference import settings
from pyreference.gene import Gene
from pyreference.mirna import MiRNA
from pyreference.pyreference_config import load_params_from_config
from pyreference.transcript import Transcript
from pyreference.utils.genomics_utils import get_unique_features_from_genomic_array_of_sets_iv, fasta_to_hash, \
    HTSeqInterval_to_feature_dict, reverse_complement, format_chrom
from pysam import FastaFile  # @UnresolvedImport
import six
import sys

CDOT_VERSION_SCHEMA = (0, 2, 0)
FASTA_LOOKUP_HAS_CHR = "chr"
FASTA_LOOKUP_NO_CHR = "no_chr"
FASTA_LOOKUP_CONTIG = "contig"
FASTA_LOOKUP = {None, FASTA_LOOKUP_HAS_CHR, FASTA_LOOKUP_NO_CHR, FASTA_LOOKUP_CONTIG}


def get_schema_version(version_tuple):
    """ Return an int which increments upon breaking changes - ie anything other than patch """
    major, minor, patch = version_tuple
    return 1000 * int(major) + int(minor)


def _load_gzip_json(gz_json_file_name, use_gzip_open=True):
    decompress_in_memory = not use_gzip_open
    if not os.path.exists(gz_json_file_name):
        raise FileNotFoundError("'%s' does not exist!" % gz_json_file_name)

    if use_gzip_open:
        try:
            with gzip.open(gz_json_file_name, "rb") as f:
                json_bytes = f.read()
        except IOError as e:
            # We sometimes get [Errno 5] Input/output error using CIFS (SMB)
            print(e, file=sys.stderr)
            if e.errno == 5:
                decompress_in_memory = True

    if decompress_in_memory:
        with open(gz_json_file_name, "rb") as f:
            gzip_bytes = f.read()
            json_bytes = gzip.GzipFile(fileobj=six.StringIO(gzip_bytes)).read()

        if use_gzip_open:
            msg = "gzip.open failed, successfully fell back on in-memory decompression\n"
            msg += "Please set use_gzip_open=False in your settings to speed up load times."
            print(msg, file=sys.stderr)

    if six.PY2:
        json_str = json_bytes
    else:
        json_str = json_bytes.decode('ascii')
    data = json.loads(json_str)

    extra_message = None
    if raw_json_version := data.get(settings.CDOT_JSON_VERSION_KEY):
        json_version = get_schema_version(raw_json_version.split("."))
        version_key = settings.CDOT_JSON_VERSION_KEY
    elif old_pyreference_version := data.get("pyreference_json_version"):
        json_version = "Old pre-cot Pyreference v%d" % old_pyreference_version
        version_key = "pyreference_json_version"
        extra_message = "PyReference switched to using cdot generated files in November 2022\n"
    else:
        raise ValueError('Invalid PyReference genes_json file: %s' % gz_json_file_name)

    cdot_schema_version = get_schema_version(CDOT_VERSION_SCHEMA)
    if cdot_schema_version != json_version:
        params = {
            "pyreference_version": pyreference.__version__,
            "cdot_schema_version": cdot_schema_version,
            "version_key": version_key,
            "json_version": json_version,
            "file_name": gz_json_file_name,
            "wiki_url": "https://github.com/SACGF/pyreference/wiki/genes_json_file",
        }
        msg = "PyReference %(pyreference_version)s requires cdot genes JSON file of schema v.%(cdot_schema_version)d\n"
        msg += "Genes JSON file '%(file_name)s' has %(version_key)s: %(json_version)s.\n"
        if extra_message:
            msg += extra_message
        msg += "Please download or re-create a genes JSON file from GTF. See %(wiki_url)s"
        raise ValueError(msg % params)

    return data


class Reference(object):
    def __init__(self, build=None, config=None, load_config_file=True, **kwargs):
        """ Construct a new reference object via:
        
            build - from pyreference config file (defaults to [global] default_build from config file) 
            config - config file (defaults to ~/pyreference.cfg)
            
            OR pass in manually:

            genome_accession
            genes_json
            trna_json
            genome_sequence_fasta
            genome_sequence_lookup
            mature_mir_sequence_fasta
            
            Any passed parameters will overwrite those from the config file
            
            stranded - interval tests are stranded? (default True) 
            
            """

        # May not need to have config file if they passed in params 
        params = {"build": build}
        config_exception = None
        try:
            if load_config_file is True:
                params = load_params_from_config(build=build, config=config)
        except OSError as e:
            config_exception = e

        # Set / Overwrite with non-null kwargs
        params.update({k: v for (k, v) in kwargs.items() if v is not None})
        self._genome_accession = params.get("genome_accession")
        self._genes_json = params.get("genes_json")
        self._trna_json = params.get("trna_json")
        self._genome_sequence_fasta = params.get("genome_sequence_fasta")
        self._genome_sequence_lookup = params.get("genome_sequence_lookup")
        self._mature_mir_sequence_fasta = params.get("mature_mir_sequence_fasta")
        self.use_gzip_open = params.get("use_gzip_open", True)
        self.stranded = params.get("stranded", True)

        # Need at least this
        REQUIRED = {
            "genome_accession": self._genome_accession,
            "genes_json": self._genes_json,
        }

        for key, data in REQUIRED.items():
            if data is None:
                message = "No '" + key + "' in"
                if kwargs:
                    six.raise_from(ValueError(message + " passed kwargs"), config_exception)
                if config_exception:
                    raise config_exception
                raise ValueError(message + " config section '%s' in file '%s'" % (params['build'], params['config']))

        if self._genome_sequence_lookup not in FASTA_LOOKUP:
            valid_values = ','.join(str(s) for s in FASTA_LOOKUP)
            raise ValueError("genome_sequence_lookup='%s' must be one of %s" % (self._genome_sequence_lookup,
                                                                                valid_values))

        self.contig_to_chrom = make_ac_name_map(self._genome_accession)
        # Store this so we can ask about config later
        self.build = params["build"]
        self._args = {"build": build, "config": config}
        self._build_params = params

    @lazy
    def _genes_dict(self):
        return _load_gzip_json(self._genes_json, self.use_gzip_open)

    def get_transcript_dict(self, transcript_id):
        """ Moves 'genome_build' down into 1st level of dict as we only need 1 """
        transcripts_by_id = self._genes_dict["transcripts"]
        tdata = transcripts_by_id[transcript_id].copy()
        genome_build = tdata.pop("genome_builds")
        tdata.update(genome_build[self._genome_accession])
        exons = tdata["exons"]
        tdata[settings.START] = exons[0][0]
        tdata[settings.END] = exons[-1][1]
        contig = tdata[settings.CONTIG]
        tdata[settings.CHROM] = self.contig_to_chrom.get(contig, contig)  # Leave as is, if not in map
        return tdata

    @lazy
    def _gene_id_lookups(self):
        gene_transcripts = defaultdict(set)
        gene_version_by_biotype = defaultdict(set)  # Set from both genes/transcripts
        for transcript_id, tdata in self._genes_dict["transcripts"].items():
            if gene_version := tdata["gene_version"]:
                gene_transcripts[gene_version].add(transcript_id)
                for biotype in tdata["biotype"]:
                    gene_version_by_biotype[biotype].add(gene_version)

        gene_version_by_symbol = {}
        for gene_version, gdata in self._genes_dict["genes"].items():
            if gene_symbol := gdata.get("gene_symbol"):
                gene_version_by_symbol[gene_symbol] = gene_version
            if biotype := gdata.get("biotype"):
                gene_version_by_biotype[biotype].add(gene_version)

        return gene_transcripts, gene_version_by_symbol, gene_version_by_biotype

    @property
    def gene_transcripts(self):
        return self._gene_id_lookups[0]

    @property
    def gene_id_by_name(self):
        return self._gene_id_lookups[1]

    @property
    def gene_ids_by_biotype(self):
        return self._gene_id_lookups[2]

    @lazy
    def genes(self):
        """ dict of {"gene_id" : Gene} """

        genes_by_id = self._genes_dict["genes"]
        genes = {}
        for gene_id in genes_by_id:
            genes[gene_id] = self.get_gene_by_id(gene_id)
        return genes

    @lazy
    def transcripts(self):
        """ dict of {"transcript_id" : Transcript} """

        # Implement in terms of genes so that it shares Transcript objects (with gene set)        
        transcripts = {}
        for g in six.itervalues(self.genes):
            for t in g.transcripts:
                transcripts[t.get_id()] = t
        return transcripts

    @lazy
    def protein_coding_genes(self):
        """Return dict of {gene_name: Gene} for protein coding genes"""
        genes_dict = {}
        for gene in self.genes_by_biotype["protein_coding"]:
            genes_dict[gene.name] = gene
        return genes_dict

    @lazy
    def genes_by_biotype(self):
        """ dict of {"biotype" : array_of_genes_biotype }
            This also includes 'tRNA' (from non-standard UCSC GTF) """

        genes_by_biotype = {}
        for (biotype, gene_ids) in self.gene_ids_by_biotype.items():
            genes = []
            for gene_id in gene_ids:
                genes.append(self.get_gene_by_id(gene_id))

            genes_by_biotype[biotype] = genes

        return genes_by_biotype

    def get_gene_by_id(self, gene_id):
        genes_by_id = self._genes_dict["genes"]
        gene_dict = genes_by_id.get(gene_id)
        if gene_dict is None:
            msg = "No Gene found with ID=%s" % gene_id
            raise ValueError(msg)

        gene_dict = gene_dict.copy()
        # Add generated transcript array
        transcripts = self.gene_transcripts.get(gene_id, [])
        gene_dict["transcripts"] = transcripts
        # Retrieve gene extents from transcript
        start = sys.maxsize
        end = 0
        chrom = None
        strand = None
        for transcript_id in transcripts:
            tdata = self.get_transcript_dict(transcript_id)
            exons = tdata["exons"]
            start = min(start, exons[0][0])
            end = max(end, exons[-1][1])
            if chrom is None:
                chrom = tdata[settings.CHROM]
                strand = tdata[settings.STRAND]

        gene_dict[settings.CHROM] = chrom
        gene_dict[settings.STRAND] = strand
        gene_dict[settings.START] = start
        gene_dict[settings.END] = end
        return Gene(self, gene_id, gene_dict)

    def get_transcript_by_id(self, transcript_id):
        transcript_dict = self.get_transcript_dict(transcript_id)
        if transcript_dict is None:
            msg = "No Transcript found with ID=%s" % transcript_dict
            raise ValueError(msg)
        return Transcript(self, transcript_id, transcript_dict)

    def get_gene_by_name(self, gene_name):
        gene_id = self.gene_id_by_name.get(gene_name)
        if gene_id is None:
            msg = "No Gene found with Name=%s" % gene_name
            raise ValueError(msg)
        return self.get_gene_by_id(gene_id)

    # @deprecated(details="Use get_gene_by_id")
    def get_gene(self, gene_id):
        return self.get_gene_by_id(gene_id)

    # @deprecated(details="Use get_transcript_by_id")
    def get_transcript(self, transcript_id):
        return self.get_transcript_by_id(transcript_id)

    def __getitem__(self, gene_symbols):
        return self.get_genes_by_name(gene_symbols)

    def get_genes_by_id(self, gene_ids):
        genes_subset = []
        for gene_id in gene_ids:
            genes_subset.append(self.get_gene_by_id(gene_id))
        return genes_subset

    def get_genes_by_name(self, gene_symbols):
        genes_subset = []
        for symbol in gene_symbols:
            genes_subset.append(self.get_gene_by_name(symbol))
        return genes_subset

    @lazy
    def mirna_mature(self):
        """Returns dict of { miRNA name: mature miR RNA sequence }"""
        if not self._mature_mir_sequence_fasta:
            msg = "Asked for miRNA mature sequence but not provided with a 'mature_mir_sequence_fasta' argument"
            raise RuntimeError(msg)
        return fasta_to_hash(self._mature_mir_sequence_fasta)

    def get_mirna(self, mirna_name):
        """Returns an MiRNA class instance for the mirna_name"""
        mirna = MiRNA(mirna_name, self.mirna_mature)
        return mirna

    @lazy
    def genome(self):
        if not self._genome_sequence_fasta:
            raise IOError("Asked for sequence but no genome_sequence_fasta provided.")
        if not os.path.exists(self._genome_sequence_fasta):
            raise IOError("Genome sequence file '%s' not found" % self.reference_fasta)

        return FastaFile(self._genome_sequence_fasta)

    def get_sequence_from_iv(self, iv, upper_case=True):
        feature_dict = HTSeqInterval_to_feature_dict(iv)
        return self.get_sequence_from_feature(feature_dict, upper_case=upper_case)

    def get_fasta_lookup_for_chrom(self, chrom):
        """ Some fasta files use contigs """

        if self._genome_sequence_lookup:
            if self._genome_sequence_lookup == FASTA_LOOKUP_HAS_CHR:
                fasta_lookup = format_chrom(chrom, want_chr=True)
            elif self._genome_sequence_lookup == FASTA_LOOKUP_NO_CHR:
                fasta_lookup = format_chrom(chrom, want_chr=False)
            elif self._genome_sequence_lookup == FASTA_LOOKUP_CONTIG:
                fasta_lookup = self.chrom_to_contig[chrom]
            else:
                raise ValueError("Unknown value for _genome_sequence_lookup: %s" % self._genome_sequence_lookup)
        else:
            fasta_lookup = chrom

        return fasta_lookup

    @lazy
    def chrom_to_contig(self):
        return {chrom: contig for contig, chrom in self.contig_to_chrom.items()}

    def get_sequence_from_feature(self, feature_dict, upper_case=True):
        """Repetitive regions are sometimes represented as lower case.
            If upper_case=True, return the sequence as upper case (Default).
            If false, do not convert case, i.e retain lower case where it was present."""

        chrom = str(feature_dict[settings.CHROM])
        start = feature_dict[settings.START]
        end = feature_dict[settings.END]
        strand = str(feature_dict[settings.STRAND])
        fasta_lookup = self.get_fasta_lookup_for_chrom(chrom)
        try:
            seq = self.genome.fetch(reference=fasta_lookup,
                                    start=start,
                                    end=end)
        except KeyError:
            self._genome_sequence_lookup

            msg = "Fasta sequence '%s' did not contain '%s'. " % (self._genome_sequence_fasta, fasta_lookup)
            if fasta_lookup != chrom:
                msg += " (converted from chrom='%s')" % chrom
            valid_values = ','.join(str(s) for s in FASTA_LOOKUP)
            params = (self._genome_sequence_lookup, valid_values, ', '.join(self.genome.references[:5]))
            msg += "You can change how chromosomes are looked up in Fasta files with 'genome_sequence_lookup'. " \
                   "Current value is '%s', allowed values = '%s'. First 5 refs in genome are %s" % params
            raise KeyError(msg)

        if strand == '-':
            seq = reverse_complement(seq)

        if upper_case:
            seq = seq.upper()
        return seq

    def get_sequence_from_features(self, features):
        sequences = []
        for feature in features:
            sequences.append(self.get_sequence_from_feature(feature))

        if sequences:
            sequence = reduce(operator.add, sequences)
        else:
            sequence = ""
        return sequence

    ###########################
    # Retrieve gene/transcript info by interval
    # Done by building a HTSeq GenomicArrayOfSets

    @lazy
    def genomic_transcripts(self):
        """ GenomicArrayOfSets containing transcripts - used for lookups """
        transcripts_gas = HTSeq.GenomicArrayOfSets("auto", self.stranded)
        for t in self.transcripts.values():
            transcripts_gas[t.iv] += t

        return transcripts_gas

    def get_transcripts_in_iv(self, iv):
        """Returns: list of transcripts in genomic interval"""
        transcripts = get_unique_features_from_genomic_array_of_sets_iv(self.genomic_transcripts, iv)
        return list(transcripts)

    def get_transcript_ids(self, iv):
        return [feature.get_id() for feature in self.get_transcripts_in_iv(iv)]

    def get_gene_names_array(self, iv):
        return list(set([t.get_gene_id() for t in self.get_transcripts_in_iv(iv)]))

    def get_gene_names(self, interval):
        """Returns a string of gene names"""
        gene_names = self.get_gene_names_array(interval)
        return " ".join(gene_names)

    ## Regions
    @lazy
    def regions(self):
        regions = HTSeq.GenomicArray("auto", stranded=self.stranded, typecode='O')

        # GenomicArray can only store 1 region for an interval, so go through 
        # and overwrite    
        for t in six.itervalues(self.transcripts):
            regions[t.iv] = "intron"

        for t in six.itervalues(self.transcripts):
            if not t.is_coding:
                for feature in t.get_features("exon"):
                    regions[feature.iv] = "non coding"

        for t in six.itervalues(self.transcripts):
            if t.is_coding:
                for utr in ["5PUTR", "3PUTR"]:
                    for feature in t.get_features(utr):
                        regions[feature.iv] = utr

        for t in six.itervalues(self.transcripts):
            if t.is_coding:
                for feature in t.get_features("CDS"):
                    regions[feature.iv] = "coding"

        return regions

    def get_regions_array(self, iv):
        """Returns: list of region types (str) in the interval"""
        _ = self.transcripts
        regions = []
        for (_, region_name) in self.regions[iv].steps():
            if region_name:
                regions.append(region_name)
        return regions

    def get_region_names(self, iv):
        region_names = set(self.get_regions_array(iv))
        return " ".join(region_names)

    def get_best_region(self, iv):
        """ Returns "best" region - if multiple pick according to order of BEST_REGION_TYPE_ORDER """

        region_names = set(self.get_regions_array(iv))
        region = None
        for r in settings.BEST_REGION_TYPE_ORDER:
            if r in region_names:
                region = r
                break
        return region

    @deprecated(details="Use get_best_region()")
    def get_region(self, iv):
        return self.get_best_region(iv)

    @lazy
    def has_chr(self):
        transcripts_by_id = self._genes_dict["transcripts"]
        some_transcript_id = six.next(six.iterkeys(transcripts_by_id))
        some_transcript = self.get_transcript_dict(some_transcript_id)
        chrom = some_transcript[settings.CHROM]
        return chrom.startswith("chr")

    def __repr__(self):
        return "PyReference (%s)" % self.build

    @lazy
    def config(self):
        params = {"build": self.build,
                  "args": self._args,
                  "reference_gtf": self._genes_dict["reference_gtf"],
                  "build_params": self._build_params.copy()}
        return params
