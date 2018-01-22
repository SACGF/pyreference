'''
Created on 22Jan.,2018


@author: dlawrence
'''
import HTSeq
import abc
from lazy import lazy

from pyreference.utils.genomics_utils import iv_from_pos_range


class GenomicRegion(object):
    ''' Base class for both Gene and Transcript '''
    def __init__(self, reference, accession_id, data_dict):
        self.reference = reference
        self.accession_id = accession_id
        self._dict = data_dict
    
    @lazy
    def iv(self):
        chrom = self._dict["chrom"]
        start = self._dict["start"]
        end = self._dict["end"]
        strand = self._dict["strand"]
        return HTSeq.GenomicInterval(chrom, start, end, strand)

    @lazy    
    def tss(self):
        return self.iv.start_d_as_pos

    def get_promoter_iv(self, promoter_range=1000):
        return iv_from_pos_range(self.tss, promoter_range)

    def get_promoter_sequence(self, promoter_range=1000):
        iv = self.get_promoter_iv(promoter_range)
        return self.reference.get_sequence_from_iv(iv)


    @abc.abstractmethod
    def get_representative_transcript(self):
        pass
