'''
Created on 19Jan.,2018

@author: dlawrence
'''
'''
TODO:

* Reference(build=XX)


'''

import lazy
import os


class Reference(object):
    def __init__(self,
                 genes_json=None,
                 trna_json=None,
                 mature_mir_fasta=None,
                 reference_fasta=None):
        
        self._genes_json = genes_json
        self._trna_json = trna_json
        self._mature_mir_fasta = mature_mir_fasta 
        self._reference_fasta = reference_fasta


    @classmethod
    def load(cls, build=None, pyreference_cfg=None):
        ''' returns a new Reference object
        
            build - from pyreference config file (default if not specified) 
            pyreference_cfg - ~/pyreference.cfg '''
        
        if pyreference_cfg is None:
            pyreference_cfg = os.path.expanduser("~/pyreference.cfg")
        
        # TODO: Verify that cfg is ok
        
        
        if build is None:
            # TODO: Get default from config
            pass
        
        kwargs = {'genes_json' : None,
                  'trna_json' : None,
                  'mature_mir_fasta' : None,
                  'reference_fasta' : None,}
        return cls(**kwargs)


    @lazy
    def genes(self):
        ''' dict of {"gene_id" : Gene} '''
        return {}

    @lazy
    def transcripts(self):
        ''' dict of {"transcript_id" : Transcript} '''
        return {}

    @lazy
    def regions(self):
        return {}
    
    @lazy
    def genes_by_biotype(self):
        return {}
    
    
    @lazy
    def genome(self):
        if not self.reference_fasta:
            raise Exception("asked for genome sequence file but Reference was not given a .fasta file")

        # @see: https://pypi.python.org/pypi/pyfasta
        key_fn = lambda key : key.split()[0] # Use first value before whitespace as keys
        return Fasta(self.reference_genome, key_fn=key_fn)
