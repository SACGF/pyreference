'''
Created on 22Jan.,2018

@author: dlawrence
'''
import HTSeq


def HTSeqInterval_to_pyfasta_feature(iv):
    return {'chr' : iv.chrom, 'start' : iv.start, 'stop' : iv.end, 'strand' : iv.strand}



def iv_from_pos_range(g_pos, range_length):
    '''p = Genomic position. 
    Returns iv 'range_length' bp upstream and 'range_length' downstream of position p'''
    return HTSeq.GenomicInterval( g_pos.chrom, g_pos.pos - range_length, g_pos.pos + range_length, g_pos.strand)
