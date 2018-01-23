from __future__ import print_function, absolute_import

from deprecated import deprecated
import gzip
import json
from lazy import lazy
import operator
import os
from pyfasta.fasta import Fasta
from six import raise_from

from pyreference.gene import Gene
from pyreference.transcript import Transcript
from pyreference.utils.genomics_utils import HTSeqInterval_to_pyfasta_feature

from .pyreference_config import load_params_from_config


def _load_gzip_json(gz_json_file_name):
    with gzip.open(gz_json_file_name) as f:
        json_bytes = f.read()
        json_str = json_bytes.decode('utf-8') 
        data = json.loads(json_str)
    return data



class Reference(object):
    def __init__(self, build=None,
                 pyreference_cfg=None,
                 **kwargs):
        ''' Construct a new reference object via:
        
            build - from pyreference config file (default if not specified) 
            pyreference_cfg - config file, ~/pyreference.cfg if not specified
            
            OR pass in the file names:
            
            genes_json
            trna_json
            mature_mir_sequence_fasta
            genome_sequence_fasta
            
            Any passed parameters will overwrite those from the config file 
            
            '''

        # May not need to have config file if they passed in params 
        config_exception = None
        try:
            params = load_params_from_config(build=build, pyreference_cfg=pyreference_cfg)
        except Exception as e:
            config_exception = e
            params = {}

        # Set / Overwrite with kwargs
        params.update(kwargs)
        self._genes_json = params.get("genes_json")
        self._trna_json = params.get("trna_json")
        self._mature_mir_sequence_fasta = params.get("mature_mir_sequence_fasta") 
        self._genome_sequence_fasta = params.get("genome_sequence_fasta")

        # Need at least this
        if self._genes_json is None:
            if kwargs:
                raise_from(ValueError("No 'genes_json' in passed kwargs"), config_exception)
                
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
        
        
        
        return {}

    @lazy
    def transcripts(self):
        ''' dict of {"transcript_id" : Transcript} '''
        
        transcripts_by_id = self._genes_dict["transcripts_by_id"]
        return {transcript_id : self.get_transcript_by_id(transcript_id) for transcript_id in transcripts_by_id}     

    @lazy
    def protein_coding_genes(self):
        '''Return dict of {gene_name: Gene} for protein coding genes'''
        genes_dict = {}
        for gene in self.genes_by_biotype["protein_coding"]:
            genes_dict[gene.name] = gene
        return genes_dict
    
#    @lazy
#    def protein_coding_genes(self):
#        ''' dict of {"gene_id" : Gene} '''
#        
#        gene_ids_by_biotype = self._genes_dict["gene_ids_by_biotype"]
#        protein_coding_gene_ids = gene_ids_by_biotype["protein_coding"] 
#        
#        return {gene_id : self.get_gene_by_id(gene_id) for gene_id in protein_coding_gene_ids}      
        

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


    @lazy
    def regions(self):
        return {}
    
    def get_gene_by_id(self, gene_id):
        genes_by_id = self._genes_dict["genes_by_id"]
        gene_dict = genes_by_id.get(gene_id)
        return Gene(self, gene_id, gene_dict)     

    def get_transcript_by_id(self, transcript_id):
        transcript_dict = self.get_transcript_dict(transcript_id)
        return Transcript(self, transcript_id, transcript_dict)


    def get_gene_by_name(self, gene_name):
        gene_id_by_name = self._genes_dict["gene_id_by_name"]
        gene_id = gene_id_by_name[gene_name]
        return self.get_gene_by_id(gene_id)
    
    
    @deprecated(reason="Use gene_gene_by_id")
    def get_gene(self, gene_id):
        return self.get_gene_by_id(gene_id)

    @deprecated(reason="Use get_transcript_by_id")
    def get_transcript(self, transcript_id):
        return self.get_transcript_by_id(transcript_id)

    
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
        ''' features: iterable of pyfasta_features '''
        
        sequences = []
        for feature in features:
            sequences.append(self.get_sequence_from_pyfasta_feature(feature))

        if sequences:
            sequence = reduce(operator.add, sequences)
        else:
            sequence = ""    
        return sequence

        




