from __future__ import print_function, absolute_import

import HTSeq
from collections import defaultdict
from lazy import lazy

from pyreference.genomic_region import GenomicRegion
from pyreference.settings import START, END, CHROM, STRAND
from pyreference.utils.genomics_utils import GenomicInterval_from_directional, dict_to_iv


class NotOnTranscriptException(Exception):
    """ For when a transcription position is not on a transcript """
    pass


class Transcript(GenomicRegion):
    def __init__(self, *args, **kwargs):
        gene = kwargs.pop("gene", None)
        super(Transcript, self).__init__(*args, **kwargs)
        self.gene = gene
    
    def get_gene_id(self):
        return self.gene.get_id() 

    @lazy
    def is_coding(self):
        return "start_codon" in self._dict

    @property
    def is_forward_strand(self):
        return self._dict["strand"] == "+"

    def get_representative_transcript(self):
        return self

    def get_features_length(self, feature_type):
        length = 0
        for feature in self.get_features_in_stranded_order(feature_type):
            length += feature[END] - feature[START] 
        return length

    def get_features(self, feature_type):
        """ returns list of HTSeq.GenomicFeature """
        genomic_features = []
        for f in self.get_features_in_stranded_order(feature_type):
            iv = dict_to_iv(f)
            genomic_feature = HTSeq.GenomicFeature(self.get_id(), feature_type, iv)
            genomic_features.append(genomic_feature)
        
        return genomic_features 

    @lazy
    def features_by_type(self):
        """ These are redundant so we re-generate them from JSON """
        fbt = defaultdict(list)

        # All in genomic order
        (left_utr, right_utr) = ("5PUTR", "3PUTR")
        if not self.is_forward_strand:  # Switch
            (left_utr, right_utr) = (right_utr, left_utr)

        cds_start = self._dict.get("cds_start")
        cds_end = self._dict.get("cds_end")

        if self.is_coding:
            left_codon_feature = {START: cds_start, END: cds_start+3}
            right_codon_feature = {START: cds_end - 3, END: cds_end}
            # cds_start/cds_end INCLUDE the start/stop codons, while the "CDS" features only includes start_codon
            cds_feature_start = cds_start
            cds_feature_end = cds_end
            if self.is_forward_strand:
                fbt["start_codon"].append(left_codon_feature)
                fbt["stop_codon"].append(right_codon_feature)
                cds_feature_end -= 3
            else:
                fbt["start_codon"].append(right_codon_feature)
                fbt["stop_codon"].append(left_codon_feature)
                cds_feature_start += 3
        else:
            cds_feature_start = None
            cds_feature_end = None

        for exon in self._dict["exons"]:  # exons in genomic order
            exon_start = exon[0]
            exon_end = exon[1]
            exon_feature = {
                START: exon_start,
                END: exon_end,
            }
            fbt["exon"].append(exon_feature)

            if self.is_coding:
                if exon_start <= cds_feature_end and exon_end >= cds_feature_start:
                    start_coding = max(cds_feature_start, exon_start)
                    stop_coding = min(cds_feature_end, exon_end)

                    cds_feature = {START: start_coding,
                                   END: stop_coding}
                    fbt["CDS"].append(cds_feature)

                if exon_start < cds_start:
                    end_non_coding = min(cds_start, exon_end)
                    utr_feature = {START: exon_start,
                                   END: end_non_coding}
                    fbt[left_utr].append(utr_feature)

                if exon_end > cds_end:
                    start_non_coding = max(cds_end, exon_start)
                    utr_feature = {START: start_non_coding,
                                   END: exon_end}
                    fbt[right_utr].append(utr_feature)

        return fbt

    def get_features_in_stranded_order(self, feature_type):
        """features returned sorted 5' -> 3' """

        is_reversed = self._dict["strand"] == '-'
        features = self.features_by_type.get(feature_type, [])
        if features:
            # Need to add this as not in there by default
            transcript_chrom = self._dict[CHROM]
            transcript_strand = self._dict[STRAND]

            for f in features:
                f[CHROM] = transcript_chrom
                f[STRAND] = transcript_strand

            features = sorted(features, key=lambda x: x[START], reverse=is_reversed)
        return features

    @lazy
    def length(self):
        return sum([exon[1] - exon[0] for exon in self._dict["exons"]])

    def get_transcript_length(self):
        return self.length

    def get_sequence_from_features(self, feature_type):
        features = self.get_features_in_stranded_order(feature_type)
        return self.reference.get_sequence_from_features(features)

    def get_transcript_sequence(self):
        return self.get_sequence_from_features("exon")

    @lazy
    def exons(self):
        """ Returns list of exon features """
        return self.get_features("exon")
    
    @lazy
    def threeputr(self):
        """ Returns the exon regions which contain 3'UTR as list of features """
        return self.get_features("3PUTR") 
    
    @lazy
    def fiveputr(self):
        """ Returns the exon regions which contain 5'UTR as list of features """
        return self.get_features("5PUTR")

    def get_coding_sequence(self):
        """ Warning: There are frame shift issues not handled here.
            Do not naively turn this into a protein - better to use existing databases """
        return self.get_sequence_from_features("CDS")

    def get_5putr_sequence(self):
        return self.get_sequence_from_features("5PUTR")

    def get_3putr_sequence(self):
        return self.get_sequence_from_features("3PUTR")

    def get_intron_ivs(self):
        """
        Purpose: Get list of intron intervals for the the transcript
        Output: A list of HTSeq Genomic Interval instances for all introns in the transcript (ordered according to strand)
        """
        intron_ivs = []
        previous_exon = None
        for exon in self.get_features("exon"):  # This is in stranded order
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
        """
        Get list of intron sequences for transcript
        Output: List of sequences (intron order and sequences are 5' to 3')
        """
        list_of_intron_intervals = self.get_intron_ivs()
        intron_sequences = []
        for intron in list_of_intron_intervals:
            intron_sequences.append(self.reference.get_sequence_from_iv(intron))
        return intron_sequences
    
    def get_genomic_position(self, pos_on_transcript):
        """
        Converts 0-based position on a transcript into 0-based position on the chromosome
        Arguments -- position relative to 5' end of transcript (int) 
        Returns -- position on chromosome (int)
        """
        #logging.debug("Searching %s for %d", self.get_id(), pos_on_transcript)
        
        running_exon_length = 0
        previous_running_exon_length = 0
        #go through exons in order, adding to running_exon_length, until you find the one with the miR seed start position in it
        for exon in self.get_features("exon"):
            running_exon_length += exon.iv.length
            if pos_on_transcript < running_exon_length: #match start is in this exon
                #logging.debug("found the exon, %r running exon length is: %r", exon, running_exon_length)
                #logging.debug("previous_running_exon_length is %r:", previous_running_exon_length)
                #logging.debug("exon start_d and end_d are: %r, %r", exon.iv.start_d, exon.iv.end_d)
                genomic_pos_of_exon_start = exon.iv.start_d
                pos_on_this_exon = pos_on_transcript - previous_running_exon_length
                #logging.debug("position on this exon is %r", pos_on_this_exon)
                if exon.iv.strand == '+':
                    genomic_position_of_match_start = genomic_pos_of_exon_start + pos_on_this_exon
                elif exon.iv.strand == '-':
                    genomic_position_of_match_start = genomic_pos_of_exon_start - pos_on_this_exon
                else:
                    raise ValueError("strand must be + or -, not: %r" % exon.iv.strand)
                return genomic_position_of_match_start
            previous_running_exon_length += exon.iv.length

        raise NotOnTranscriptException("%s didn't contain %s" % (self.get_id(), pos_on_transcript))

    def __repr__(self):
        coding_str = " (coding)" if self.is_coding else ""
        return "Transcript %s%s: length %d" % (self.get_id(), coding_str, self.length)
