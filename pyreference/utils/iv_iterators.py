"""
Created on 22Jan.,2018

@author: dlawrence
"""
import HTSeq
import os

from pyreference.utils.genomics_utils import format_chrom


def format_chrom_iterator(iterator, want_chr):
    for iv in iterator:
        iv.chrom = format_chrom(iv.chrom, want_chr)
        yield iv


def load_iv_iterator(file_name):
    suffix = os.path.splitext(file_name)[1].lower()
    iterator = None
    if suffix == ".bam":
        iterator = bam_iv_iterator(file_name)
    elif suffix == ".bed":
        # TODO:
        raise NotImplementedError("No bed reader yet")
    elif suffix == '.sam':
        iterator = sam_iv_iterator(file_name)
    elif suffix in ['.gff', '.gtf']:
        iterator = gff_iv_iterator(file_name)
    else:
        raise ValueError("Unknown input_file_type of " + suffix)
    return iterator


def chromosome_filter_iterator(chromosomes, iterator):
    for iv in iterator:
        if iv.chrom in chromosomes:
            yield iv


def bam_iv_iterator(bam_file):
    for aln in HTSeq.BAM_Reader(bam_file):
        if aln.aligned:
            yield aln.iv


def sam_iv_iterator(sam_file):
    for aln in HTSeq.SAM_Reader(sam_file):
        if aln.aligned:
            yield aln.iv


def gff_iv_iterator(gtf_file):
    for feature in HTSeq.GFF_Reader(gtf_file):
        yield feature.iv
