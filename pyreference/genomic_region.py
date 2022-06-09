from __future__ import print_function, absolute_import

import abc
import six
from lazy import lazy

from pyreference.utils.genomics_utils import iv_from_pos_range, \
    iv_from_pos_directional_before_after, dict_to_iv


class GenomicRegion(object):
    """ Base class for both Gene and Transcript """
    def __init__(self, reference, accession_id, data_dict):
        self.reference = reference
        self.accession_id = accession_id
        self._dict = data_dict
    
    def get_id(self):
        return self.accession_id

    @property
    def biotype(self):
        return '/'.join(sorted(self.get_biotypes()))

    def get_biotypes(self):
        # On gene it's a string
        biotype = self._dict["biotype"]
        if isinstance(biotype, six.string_types):
            biotypes = biotype.split(",")
        elif isinstance(biotype, list):
            biotypes = biotype
        return biotypes
    
    @lazy
    def iv(self):
        return dict_to_iv(self._dict)

    @lazy    
    def tss(self):
        """ (Representative) Transcript Start Site
            This is NOT the most 5' position (use iv.start_d_as_pos for that) """
        transcript_iv = self.get_representative_transcript().iv
        return transcript_iv.start_d_as_pos

    def get_promoter_iv(self, promoter_range=1000):
        return iv_from_pos_range(self.tss, promoter_range)

    def get_promoter_sequence(self, promoter_range=1000):
        iv = self.get_promoter_iv(promoter_range)
        return self.reference.get_sequence_from_iv(iv)

    def get_promoter_iv_custom_range(self, upstream_distance, downstream_distance):
        """Get any interval surrounding TSS
        Note: total length of interval = upstream_distance + downstream_distance (The TSS base is included in downstream_distance)"""
        return iv_from_pos_directional_before_after(self.tss, upstream_distance, downstream_distance)
    
    def get_promoter_sequence_custom_range(self, upstream_distance, downstream_distance):
        iv = self.get_promoter_iv_custom_range(upstream_distance, downstream_distance)
        return self.reference.get_sequence_from_iv(iv)

    @abc.abstractmethod
    def get_representative_transcript(self):
        pass
