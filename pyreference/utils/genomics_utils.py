"""
Created on 22Jan.,2018

@author: dlawrence
"""

from Bio.Seq import Seq
from Bio import SeqIO
import HTSeq

try:
    from sys import intern
except (ImportError,AttributeError):
    pass
    
from pyreference.settings import CHROM, START, END, STRAND


def HTSeqInterval_to_feature_dict(iv):
    return {CHROM: iv.chrom, START: iv.start, END: iv.end, STRAND: iv.strand}


def dict_to_iv(data):
    chrom = str(data[CHROM])
    start = data[START]
    end = data[END]
    strand = str(data[STRAND])
    return HTSeq.GenomicInterval(chrom, start, end, strand)


def iv_from_pos_range(g_pos, range_length):
    """p = Genomic position.
    Returns iv 'range_length' bp upstream and 'range_length' downstream of position p"""
    return HTSeq.GenomicInterval( g_pos.chrom, g_pos.pos - range_length, g_pos.pos + range_length, g_pos.strand)


def iv_from_pos_directional_before_after(g_pos, upstream_length, downstream_length):
    """Note: The g_pos base is assumed to be included in downstream_length
    e.g upstream_length=100, downstream_length=100 has total length=200 bp
    Check that this function does what you're expecting - setting upstream length to 0 will change the position of start_d on neg strand..."""
    if g_pos.strand == '+':
        start = g_pos.pos - upstream_length
        end = g_pos.pos + downstream_length
    elif g_pos.strand == '-':
        start = g_pos.pos - downstream_length
        end = g_pos.pos + upstream_length
    else:
        raise ValueError("Unknown strand in genomic position %s" % g_pos)

    return HTSeq.GenomicInterval( g_pos.chrom, start, end, g_pos.strand)


def GenomicInterval_from_directional( chrom, start_d, length, strand="." ):
    """ Fix bug in HTSeq:
        HTSeq.GenomicInterval_from_directional throws 'str' object has no attribute 'se' """
    strand = intern( strand )
    if strand != "-":
        return HTSeq.GenomicInterval( chrom, start_d, start_d+length, strand )
    else:
        return HTSeq.GenomicInterval( chrom, start_d-length+1, start_d+1, strand )


def get_unique_features_from_genomic_array_of_sets_iv(genomic_array_of_sets, iv):
    """Collapse genomic array of sets into a unique list of the values in that region"""
    all_values = set()
    for iv, values_in_iv in genomic_array_of_sets[iv].steps():
        all_values.update(values_in_iv)
    return all_values


# HTSeq documentation says end = the end of the interval. Following Python convention for ranges, this
# is *ONE MORE* than the coordinate of the last base that is considered part of the sequence.
def last_base(iv):
    last_base_pos = iv.end_as_pos
    last_base_pos.pos -= 1
    return last_base_pos


def opposite_strand(strand):
    opposites = {"+": "-",
                 "-": "+"}
    o = opposites.get(strand)
    if o is None:
        raise ValueError("Unknown strand '%s'" % strand)
    return o


def format_chrom(chrom, want_chr):
    """ Pass in a chromosome (unknown format), return in your format
        @param chrom: chromosome ID (eg 1 or 'chr1')
        @param want_chr: Boolean - whether you want "chr" at the beginning of chrom
        @return: "chr1" or "1" (for want_chr True/False) 
    """
    if want_chr:   
        if chrom.startswith("chr"):
            return chrom
        else:
            return "chr%s" % chrom
    else:
        return chrom.replace("chr", "")


def fasta_to_hash(fasta):
    indexed_fasta = {}
    for record in SeqIO.parse(fasta, "fasta"):
        indexed_fasta[record.id] = str(record.seq)
    return indexed_fasta


def reverse_complement(dna_sequence):
    seq = Seq(dna_sequence)
    return str(seq.reverse_complement())
