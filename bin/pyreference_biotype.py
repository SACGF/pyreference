#!/usr/bin/env python
"""
Created on 22Jan.,2018

@author: dlawrence
"""

from __future__ import print_function, absolute_import

import HTSeq
import abc
from argparse import ArgumentParser
from collections import Counter, defaultdict
import csv
import logging
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.cm import ScalarMappable
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import os
from pyreference import Reference
from pyreference.utils import iv_iterators
from pyreference.utils.file_utils import name_from_file_name, file_or_file_name, \
    mk_path_for_file
from pyreference.utils.genomics_utils import opposite_strand, format_chrom
import six
import seaborn as sns
import pandas as pd
import numpy as np

MAX_READ_LENGTH = 200


def handle_args():
    parser = ArgumentParser(description='Collect stats on read length and biotype')
    parser.add_argument("--intervals", help='.bed/.gtf etc file')
    parser.add_argument("--intervals-name", help="Used in graphs")
    parser.add_argument("--reverse-strand", action='store_true', help="Reverse strand before testing region")
    parser.add_argument("bam")
    return parser.parse_args()


def get_length_counts(bam, regions_array, has_chr, reverse_strand):
    length_counters = defaultdict(Counter)

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


def create_biotype_regions_array(reference, interesting_biotypes=None):
    """ genes_by_biotype : dict of {"biotype" : genes[]}
        interesting_biotypes : List of Strings corresponding to biotype keys (everything else is 'other') """

    # In HTSeq v1.99.2 "auto" GenomicArrays create non-infinite chromosome arrays if 1st accessed via a set
    # so you can get IndexError: stop too large accessing the array later, see https://github.com/htseq/htseq/issues/38
    chromosomes = set()

    if interesting_biotypes is None:
        interesting_biotypes = ['protein_coding', 'rRNA', 'lincRNA', 'misc_RNA', 'snRNA', 'miRNA', 'snoRNA', 'tRNA']

    other_biotypes = []  # Everything non interesting is counted as "other"
    for biotype in reference.genes_by_biotype.keys():
        if biotype not in interesting_biotypes:
            other_biotypes.append(biotype)

    def get_biotype(gene):
        if gene.biotype in other_biotypes:
            return "other"
        elif gene.biotype == "misc_RNA":
            if gene.name:
                if gene.name.startswith("RNY"):
                    return "yRNA"
        return gene.biotype

    regions = HTSeq.GenomicArray("auto", stranded=True, typecode='O')
    for transcript in six.itervalues(reference.transcripts):
        # Antisense: Read is in the region of a transcript, but on the opposite strand.
        antisense_iv = transcript.iv.copy()
        antisense_iv.strand = opposite_strand(antisense_iv.strand)

        # This should make all chroms as we're iterating through all transcripts above
        if antisense_iv.chrom not in regions.chrom_vectors:
            regions.add_chrom(antisense_iv.chrom)
        regions[antisense_iv] = "anti-sense"

    for gene in six.itervalues(reference.genes):
        for transcript in gene.transcripts:
            regions[transcript.iv] = "introns"

    for biotype in other_biotypes + interesting_biotypes:
        genes = reference.genes_by_biotype.get(biotype, [])
        for gene in genes:
            for t in gene.transcripts:
                for exon in t.get_features("exon"):
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
    
    #Set up  
    sample_name = name_from_file_name(args.bam)
    csv_file = "%s.read_counts.regions.csv" % sample_name
    graph_image =  "%s.read_counts.regions.png" % sample_name
    
    # Use this as a test platform to load reference
    reference = Reference()
    print("Reference is", reference) #To confirm the annotation you're using is what you intended. 
    
    regions_array = create_biotype_regions_array(reference)

    if args.intervals:
        add_intervals_to_regions_array(regions_array, args.intervals, args.intervals_name, reference.has_chr,
                                       args.reverse_strand)

    print("Loaded reference - reading BAM")
    length_counters = get_length_counts(args.bam, regions_array, reference.has_chr, args.reverse_strand)

    # Find min/max position used
    start = six.MAXSIZE
    end = 0
    for i, c in length_counters.items():
        if len(c):
            start = min(start, i)
            end = max(end, i)

    biotype_colors = {'other': "red",
                      'protein_coding': "DarkGreen",
                      'rRNA': "orange",
                      'lincRNA': "purple",
                      'misc_RNA': "brown",
                      'snRNA': "yellow",
                      'miRNA': "blue",
                      'snoRNA': "cyan",
                      "unaligned": "black",
                      "tRNA": "magenta",
                      "intergenic": "grey",
                      "introns": "pink",
                      "anti-sense": "silver",
                      "yRNA": "Coral"}
    if args.intervals:
        biotype_colors[args.intervals_name] = "lightgreen"

    largs = range(start, end + 1)
    labels = sorted(biotype_colors.keys())
    colors = []
    arrays = []
    total_counts = Counter()

    for k in labels:
        colors.append(biotype_colors[k])
        counts_for_lengths = np.zeros(len(largs), dtype='i')
        for l in largs:
            length_count = length_counters[l][k]
            counts_for_lengths[l - start] = length_count
            total_counts[k] += length_count
        arrays.append(counts_for_lengths)
        
    df = pd.DataFrame(data=arrays, columns=largs, index=labels)
    df = df.T
    df.to_csv(csv_file)

    ### Graph data ###
    sns.set_theme(context='paper', style="ticks", font_scale=1.1)
    legend_kwargs = {'loc' : 'center left',
                     'prop' : {'size': 8.5},
                     'bbox_to_anchor' : (1.01, 0.5)}
    
#     df = pd.read_csv(csv_file, index_col=0)
    print(df)
     
    fig = Figure(dpi=300, figsize=(4.8, 3.1))
    
    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(111)
     
    #Make stacked bar chart
    bottom = np.zeros(len(df.index), dtype='i')
    for i, label in enumerate(df.columns):
        counts = df[label]
        color = biotype_colors[label]
        _ = ax.bar(df.index, counts, label=label, color=color, bottom=bottom, linewidth=0)
        bottom += counts
    
    #Format chart
    ax.set_xlabel("Length (nt)")
    ax.set_ylabel("Read counts")
     
    ax.set_ylim(ymin=0)
    ax.set_xlim(xmin=(min(largs) - 0.7), xmax=(max(largs) + 0.7))
    
    #Shrink to fit legend
    fig.tight_layout(rect=[0,0,0.7,1]) #left, bottom, right, top
     
    ax.legend(**legend_kwargs)
    mk_path_for_file(graph_image)
     
    canvas = FigureCanvasAgg(fig)
    canvas.print_png(graph_image)

if __name__ == '__main__':
    main()
