from __future__ import print_function, absolute_import

import HTSeq
from deprecated import deprecated
import gzip
import json
from lazy import lazy
import operator
import os
from pyfasta.fasta import Fasta
import six

from pyreference import settings
from pyreference.gene import Gene
from pyreference.mirna import MiRNA
from pyreference.pyreference_config import load_params_from_config
from pyreference.settings import BEST_REGION_TYPE_ORDER
from pyreference.transcript import Transcript
from pyreference.utils.genomics_utils import HTSeqInterval_to_pyfasta_feature, \
    get_unique_features_from_genomic_array_of_sets_iv, fasta_to_hash


def _load_gzip_json(gz_json_file_name):
    with gzip.open(gz_json_file_name) as f:
        json_bytes = f.read()
        if six.PY2:
            json_str = json_bytes
        else:
            json_str = json_bytes.decode('ascii')
        data = json.loads(json_str)
        
    pyreference_json_version = data[settings.PYREFERENCE_JSON_VERSION_KEY]
    if settings.PYREFERENCE_JSON_VERSION != pyreference_json_version:
        params = {"version_key" : settings.PYREFERENCE_JSON_VERSION_KEY,
                  "current_version" : settings.PYREFERENCE_JSON_VERSION,
                  "json_version" : pyreference_json_version,
                  "file_name" : gz_json_file_name}
        msg = "PyReference with %(version_key)s %(current_version)d attempted to load '%(file_name)s' with %(version_key)s: %(json_version)d.\n" % params
        msg +=  "Please re-create with this version of pyreference_gtf_to_json.py."
        raise ValueError(msg)
      
    return data


class Reference(object):
    def __init__(self, build=None, config=None, **kwargs):
        ''' Construct a new reference object via:
        
            build - from pyreference config file (defaults to [global] default_build from config file) 
            config - config file (defaults to ~/pyreference.cfg)
            
            OR pass in the file names:
            
            genes_json
            trna_json
            genome_sequence_fasta
            mature_mir_sequence_fasta
            
            
            Any passed parameters will overwrite those from the config file 
            
            stranded - interval tests are stranded? (default True) 
            
            '''

        # May not need to have config file if they passed in params 
        config_exception = None
        try:
            params = load_params_from_config(build=build, config=config)
        except Exception as e:
            config_exception = e
            params = {}

        # Set / Overwrite with non-null kwargs
        params.update({k : v for (k,v) in kwargs.items() if v is not None})
        self._genes_json = params.get("genes_json")
        self._trna_json = params.get("trna_json")
        self._genome_sequence_fasta = params.get("genome_sequence_fasta")
        self._mature_mir_sequence_fasta = params.get("mature_mir_sequence_fasta") 
        self.stranded = kwargs.get("stranded", True)

        # Need at least this
        if self._genes_json is None:
            if kwargs:
                six.raise_from(ValueError("No 'genes_json' in passed kwargs"), config_exception)
            raise config_exception


    @lazy
    def _genes_dict(self):
        return _load_gzip_json(self._genes_json)

    def get_transcript_dict(self, transcript_id):
        transcripts_by_id = self._genes_dict["transcripts_by_id"]
        return transcripts_by_id[transcript_id]


    @lazy
    def genes(self):
        ''' dict of {"gene_id" : Gene} '''

        genes_by_id = self._genes_dict["genes_by_id"]
        genes = {}
        for gene_id in genes_by_id:
            genes[gene_id] = self.get_gene_by_id(gene_id)
        return genes
    
    @lazy
    def transcripts(self):
        ''' dict of {"transcript_id" : Transcript} '''

        # Implement in terms of genes so that it shares Transcript objects (with gene set)        
        transcripts = {}
        for g in six.itervalues(self.genes):
            for t in g.transcripts:
                transcripts[t.get_id()] = t
        return transcripts

    @lazy
    def protein_coding_genes(self):
        '''Return dict of {gene_name: Gene} for protein coding genes'''
        genes_dict = {}
        for gene in self.genes_by_biotype["protein_coding"]:
            genes_dict[gene.name] = gene
        return genes_dict

    @lazy
    def genes_by_biotype(self):
        ''' dict of {"biotype" : array_of_genes_biotype }
            This also includes 'tRNA' (from non-standard UCSC GTF) '''
        gene_ids_by_biotype = self._genes_dict["gene_ids_by_biotype"]

        genes_by_biotype = {}
        for (biotype, gene_ids) in gene_ids_by_biotype.items():
            genes = []
            for gene_id in gene_ids:
                genes.append(self.get_gene_by_id(gene_id))
                
            genes_by_biotype[biotype] = genes
            
        return genes_by_biotype


    def get_gene_by_id(self, gene_id):
        genes_by_id = self._genes_dict["genes_by_id"]
        gene_dict = genes_by_id.get(gene_id)
        if gene_dict is None:
            msg = "No Gene found with ID=%s" % gene_id
            raise ValueError(msg)
        return Gene(self, gene_id, gene_dict)     

    def get_transcript_by_id(self, transcript_id):
        transcript_dict = self.get_transcript_dict(transcript_id)
        if transcript_dict is None:
            msg = "No Transcript found with ID=%s" % transcript_dict
            raise ValueError(msg)
        return Transcript(self, transcript_id, transcript_dict)


    def get_gene_by_name(self, gene_name):
        gene_id_by_name = self._genes_dict["gene_id_by_name"]
        gene_id = gene_id_by_name.get(gene_name)
        if gene_id is None:
            msg = "No Gene found with Name=%s" % gene_name
            raise ValueError(msg)
        return self.get_gene_by_id(gene_id)
    
    
    @deprecated(reason="Use get_gene_by_id")
    def get_gene(self, gene_id):
        return self.get_gene_by_id(gene_id)

    @deprecated(reason="Use get_transcript_by_id")
    def get_transcript(self, transcript_id):
        return self.get_transcript_by_id(transcript_id)


    @lazy
    def mirna_mature(self):
        '''Returns dict of { miRNA name: mature miR RNA sequence }'''
        if not self._mature_mir_sequence_fasta:
            msg = "Asked for miRNA mature sequence but not provided with a 'mature_mir_sequence_fasta' argument"
            raise RuntimeError(msg)
        return fasta_to_hash(self._mature_mir_sequence_fasta)

    def get_mirna(self, mirna_name):
        '''Returns an MiRNA class instance for the mirna_name'''
        mirna = MiRNA(mirna_name, self.mirna_mature)
        return mirna


    @lazy
    def genome(self):
        if not self._genome_sequence_fasta:
            raise IOError("Asked for sequence but no genome_sequence_fasta provided.")
        if not os.path.exists(self._genome_sequence_fasta):
            raise IOError("Genome sequence file '%s' not found" % self.reference_fasta)

        # @see: https://pypi.python.org/pypi/pyfasta
        key_fn = lambda key : key.split()[0] # Use first value before whitespace as keys
        return Fasta(self._genome_sequence_fasta, key_fn=key_fn)

    def get_sequence_from_iv(self, iv, upper_case=True):
        pyfasta_feature = HTSeqInterval_to_pyfasta_feature(iv)
        return self.get_sequence_from_pyfasta_feature(pyfasta_feature, upper_case=upper_case)
        

    def get_sequence_from_pyfasta_feature(self, pyfasta_feature, upper_case=True):
        '''Repetitive regions are sometimes represented as lower case.
            If upper_case=True, return the sequence as upper case (Default).
            If false, do not convert case, i.e retain lower case where it was present.'''

        seq = self.genome.sequence(pyfasta_feature, one_based=False)
        if upper_case:    
            seq = seq.upper()
        return seq


    def get_sequence_from_features(self, features):
        sequences = []
        for feature in features:
            sequences.append(self.get_sequence_from_pyfasta_feature(feature))

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
        ''' GenomicArrayOfSets containing transcripts - used for lookups '''
        transcripts_gas = HTSeq.GenomicArrayOfSets("auto", self.stranded)
        for t in self.transcripts.values():
            transcripts_gas[t.iv] += t

        return transcripts_gas

    def get_transcripts_in_iv(self, iv):
        '''Returns: list of transcripts in genomic interval'''
        transcripts = get_unique_features_from_genomic_array_of_sets_iv(self.genomic_transcripts, iv)
        return list(transcripts)

    def get_transcript_ids(self, iv):
        return [feature.name for feature in self.get_transcripts_in_iv(iv)]

    def get_gene_names_array(self, iv):
        return list(set([t.get_gene_id() for t in self.get_transcripts_in_iv(iv)]))

    def get_gene_names(self, interval):
        '''Returns a string of gene names'''
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
        '''Returns: list of region types (str) in the interval'''
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
        ''' Returns "best" region - if multiple pick according to order of BEST_REGION_TYPE_ORDER '''

        region_names = set(self.get_regions_array(iv))
        region = None
        for r in BEST_REGION_TYPE_ORDER:
            if r in region_names:
                region = r
                break
        return region

    @deprecated(reason="Use get_best_region()")
    def get_region(self, iv):
        return self.get_best_region(iv)


    @lazy
    def has_chr(self):
        transcripts_by_id = self._genes_dict["transcripts_by_id"]
        some_transcript = six.next(six.itervalues(transcripts_by_id))
        chrom = some_transcript["chr"]
        return chrom.startswith("chr")

