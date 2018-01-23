'''
Created on 19Jan.,2018

TODO: ReferenceArgs to handle changing etc


@author: dlawrence
'''

import HTSeq
from argparse import ArgumentParser
from collections import Counter
import six
import sys

import numpy as np
from pyreference import Reference, reference
from pyreference.utils import iv_iterators
from pyreference.utils.csv_utils import write_csv_dict
from pyreference.utils.file_utils import name_from_file_name
from pyreference.utils.genomics_utils import opposite_strand, format_chrom
from pyreference.utils.graphing import write_stacked_bar_graph, \
    write_individual_bar_graphs


MAX_READ_LENGTH = 200

def handle_args():
    parser = ArgumentParser(description='Collect stats on read length and biotype')
    parser.add_argument("--intervals", help='.bed/.gtf etc file')
    parser.add_argument("--intervals-name", help="Used in graphs")
    parser.add_argument("--reverse-strand", action='store_true',  help="Reverse strand before testing region")
    parser.add_argument("bam")
    return parser.parse_args()

def get_length_counts(bam, regions_array, has_chr, reverse_strand):
    length_counters = [Counter() for _ in range(MAX_READ_LENGTH + 1)]

    for aln in HTSeq.BAM_Reader(bam):
        length = len(aln.read)
        read_region = None
        if aln.aligned:
            aln.iv.chrom = format_chrom(aln.iv.chrom, has_chr)
            if reverse_strand:
                aln.iv.strand = opposite_strand(aln.iv.strand)

            region_overlap_length = 0
            for iv, r in regions_array[aln.iv].steps():
                if iv.length > region_overlap_length:
                    read_region = r
                    region_overlap_length = iv.length
            
            if read_region is None:
                read_region = "intergenic"
        else:
            read_region = "unaligned"
        length_counters[length][read_region] += 1
    
    return length_counters


def create_biotype_regions_array(genes_by_biotype, interesting_biotypes=None):
    ''' genes_by_biotype : dict of {"biotype" : genes[]}
        interesting_biotypes : List of Strings corresponding to biotype keys (everything else is 'other') '''

    if interesting_biotypes is None:
        interesting_biotypes = ['protein_coding', 'rRNA', 'lincRNA', 'misc_RNA', 'snRNA', 'miRNA','snoRNA', 'tRNA']
    
    other_biotypes = [] # Everything non interesting is counted as "other"
    for biotype in genes_by_biotype.keys():
        if biotype not in interesting_biotypes:
            other_biotypes.append(biotype)

    def get_biotype(gene):
        if gene.biotype in other_biotypes:
            return "other"
        elif gene.biotype == "misc_RNA":
            if gene.name.startswith("RNY"):
                return "yRNA"
        return gene.biotype

    regions = HTSeq.GenomicArray( "auto", stranded=True, typecode='O' )
    for transcript in six.itervalues(reference.transcripts): #@UndefinedVariable
        #Antisense: Read is in the region of a transcript, but on the opposite strand.
        antisense_iv = transcript.iv.copy()
        antisense_iv.strand = opposite_strand(antisense_iv.strand)
        regions[antisense_iv] = "anti-sense"
    
    for gene in genes_by_biotype["protein_coding"]:
        for transcript in gene.transcripts:
            regions[transcript.iv] = "introns" 

    for biotype in other_biotypes + interesting_biotypes:
        for gene in genes_by_biotype[biotype]:
            for t in gene.transcripts:
                for exon in t.features_by_type["exon"]:
                    regions[exon.iv] = get_biotype(gene)
    return regions


def add_intervals_to_regions_array(regions_array, intervals, intervals_name, has_chr, reverse_strand):
    intervals_name = intervals_name or name_from_file_name(intervals)
    
    iv_iterator = iv_iterators.load_iv_iterator(intervals)
    for iv in iv_iterators.format_chrom_iterator(iv_iterator, has_chr):
        # Fix for 'Non-stranded index used for stranded GenomicArray.'
        if iv.strand == '.':
            iv.strand = '+'
            regions_array[iv] = intervals_name
            iv.strand = '-'
            regions_array[iv] = intervals_name
        else:
            regions_array[iv] = intervals_name



def main():
    args = handle_args()

    # Use this as a test platform to load reference
    reference = Reference()

    regions_array = create_biotype_regions_array(reference.genes_by_biotype)

    if args.intervals:
        add_intervals_to_regions_array(regions_array, args.intervals, args.intervals_name, reference.has_chr, args.reverse_strand)
    
    length_counters = get_length_counts(args.bam, regions_array, reference.has_chr, args.reverse_strand)
        
    start = sys.maxint
    end = 0
    for (i, c) in enumerate(length_counters):
        if len(c):
            start = min(start, i)
            end = max(end, i)
    
    biotype_colors = {'other' : "red",
                      'protein_coding' : "DarkGreen",
                      'rRNA' : "orange",
                      'lincRNA' : "purple",
                      'misc_RNA' : "brown",
                      'snRNA' : "yellow",
                      'miRNA' : "blue",
                      'snoRNA' : "cyan",
                      "unaligned" : "black",
                      "tRNA" : "magenta",
                      "intergenic" : "grey",
                      "introns" : "pink",
                      "anti-sense" : "silver",
                      "yRNA" : "Coral",
    }
    if args.intervals:
        biotype_colors[args.intervals_name] = "lightgreen"
    
    largs = range(start, end+1)
    labels = sorted(biotype_colors.keys())
    colors = []
    arrays = []
    total_counts = Counter()
    
    for k in labels:
        colors.append(biotype_colors[k])
        counts_for_lengths = np.zeros(len(largs), dtype='i')
        for l in largs:
            length_count = length_counters[l][k]
            counts_for_lengths[l-start] = length_count
            total_counts[k] += length_count
        arrays.append(counts_for_lengths)
    
    graph_types = {"read_counts" : lambda x : x}
    
    sample_name = name_from_file_name(args.bam)
    
    for (y_label, op) in graph_types.iteritems():
        graph_image =  "%s.%s.regions.png" % (sample_name, y_label)
        data = [op(i) for i in arrays]
        
        ncol = max(len(labels), 5)
        legend_kwargs = {'loc' : 'center left',
                         'prop' : {'size': 8},
                         'bbox_to_anchor' : (1.01, 0.5)}
        write_stacked_bar_graph(graph_image, largs, data, labels, colors, y_label=y_label, legend_side_ratio=0.15, legend_kwargs=legend_kwargs)
        write_individual_bar_graphs(graph_image, largs, data, labels, colors, y_label=y_label)
    
    # Write region stats CSV
    csv_file = "%s.biotype_stats.csv" % sample_name
    write_csv_dict(csv_file, labels, [total_counts])




if __name__ == '__main__':
    main()