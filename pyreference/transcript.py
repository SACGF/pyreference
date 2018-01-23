from __future__ import print_function, absolute_import

import HTSeq
from deprecated import deprecated
from lazy import lazy
import logging

from pyreference.genomic_region import GenomicRegion
from pyreference.settings import START, END, IS_CODING, CHROM, STRAND
from pyreference.utils.genomics_utils import GenomicInterval_from_directional, dict_to_iv


class NotOnTranscriptException(Exception):
    ''' For when a transcription position is not on a transcript '''
    pass


class Transcript(GenomicRegion):
    def __init__(self, *args, **kwargs):
        gene = kwargs.pop("gene", None)
        super(Transcript, self).__init__(*args, **kwargs)
        self.gene = gene
    
    def get_gene_id(self):
        return self.gene.get_id() 

    @property
    def is_coding(self):
        return self._dict[IS_CODING]

    def get_representative_transcript(self):
        return self


    def get_features_length(self, feature_type):
        length = 0
        for feature in self.get_features_in_stranded_order(feature_type):
            length += feature[END] - feature[START] 
        return length

    #@deprecated(reason="Use get_features_in_stranded_order")
    def get_features(self, feature_type):
        genomic_features = []
        for f in self.get_features_in_stranded_order(feature_type):
            iv = dict_to_iv(f)
            genomic_feature = HTSeq.GenomicFeature(self.get_id(), feature_type, iv)
            genomic_features.append(genomic_feature)
        
        return genomic_features 


    def get_features_in_stranded_order(self, feature_type):
        '''features returned sorted 5' -> 3' '''

        is_reversed = self._dict["strand"] == '-'
        if is_reversed:
            stranded_start = END
        else:
            stranded_start = START
            
        features_by_type = self._dict["features_by_type"]
        features = features_by_type.get(feature_type, [])
        if features:
            # Need to add this as not in there by default
            transcript_chrom = self._dict[CHROM]
            transcript_strand = self._dict[STRAND]

            for f in features:
                f[CHROM] = transcript_chrom
                f[STRAND] = transcript_strand

            features = sorted(features, key=lambda x : x[stranded_start], reverse=is_reversed)
        return features

    @lazy
    def length(self):
        return self.get_features_length("exon")
    
    @deprecated(reason="Use Transcript.length")
    def get_transcript_length(self):
        return self.length

    def get_sequence_from_features(self, feature_type):
        features = self.get_features_in_stranded_order(feature_type)
        return self.reference.get_sequence_from_features(features)

    def get_transcript_sequence(self):
        return self.get_sequence_from_features("exon")
    
    def get_coding_sequence(self):
        ''' Warning: There are frame shift issues not handled here.
            Do not naively turn this into a protein - better to use existing databases '''
        return self.get_sequence_from_features("CDS")

    def get_5putr_sequence(self):
        return self.get_sequence_from_features("5PUTR")

    def get_3putr_sequence(self):
        return self.get_sequence_from_features("3PUTR")

    def get_intron_ivs(self):
        '''
        Purpose: Get list of intron intervals for the the transcript
        Output: A list of HTSeq Genomic Interval instances for all introns in the transcript (ordered according to strand)
        '''
        intron_ivs = []
        previous_exon = None
        for exon in self.get_features("exon"): # This is in stranded order
            if previous_exon:
                # HTSeq ends are 1 past the last base of the sequence.
                # Thus for touching sequences like exons/introns, first_seq.end = second_seq.start
                intron_start = previous_exon.iv.end_d
                intron_end = exon.iv.start_d
                intron_length = abs(intron_end - intron_start)
                intron = GenomicInterval_from_directional(exon.iv.chrom, intron_start, intron_length, exon.iv.strand)
                intron_ivs.append(intron)
            previous_exon = exon
        return intron_ivs
    
    def get_intron_sequences(self):
        '''
        Get list of intron sequences for transcript
        Output: List of sequences (intron order and sequences are 5' to 3')
        '''
        list_of_intron_intervals = self.get_intron_ivs()
        intron_sequences = []
        for intron in list_of_intron_intervals:
            intron_sequences.append(self.reference.get_sequence_from_iv(intron))
        return intron_sequences
    
    
    def get_genomic_position(self, pos_on_transcript):
        '''
        Converts 0-based position on a transcript into 0-based position on the chromosome
        Arguments -- position relative to 5' end of transcript (int) 
        Returns -- position on chromosome (int)
        '''    
        logging.debug("Searching %s for %d", self.get_id(), pos_on_transcript)
        
        running_exon_length = 0
        previous_running_exon_length = 0
        #go through exons in order, adding to running_exon_length, until you find the one with the miR seed start position in it
        for exon in self.get_features("exon"):
            running_exon_length += exon.iv.length
            if pos_on_transcript < running_exon_length: #match start is in this exon
                logging.debug("found the exon, %r running exon length is: %r", exon, running_exon_length)
                logging.debug("previous_running_exon_length is %r:", previous_running_exon_length)
                logging.debug("exon start_d and end_d are: %r, %r", exon.iv.start_d, exon.iv.end_d)
                genomic_pos_of_exon_start = exon.iv.start_d
                pos_on_this_exon = pos_on_transcript - previous_running_exon_length
                logging.debug("position on this exon is %r", pos_on_this_exon)
                if exon.iv.strand == '+':
                    genomic_position_of_match_start = genomic_pos_of_exon_start + pos_on_this_exon
                elif exon.iv.strand == '-':
                    genomic_position_of_match_start = genomic_pos_of_exon_start - pos_on_this_exon
                else:
                    raise ValueError("strand must be + or -, not: %r" % exon.iv.strand)
                return genomic_position_of_match_start
            previous_running_exon_length += exon.iv.length

        raise NotOnTranscriptException("%s didn't contain %s" % (self.get_id(), pos_on_transcript))

