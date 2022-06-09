#!/usr/bin/env python

from __future__ import print_function, absolute_import
from argparse import ArgumentParser
from collections import Counter, defaultdict
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from pyreference import Reference, ReferenceArgumentParser
from pyreference.utils import iv_iterators
from pyreference.utils.file_utils import name_from_file_name, mk_path_for_file
from pyreference.utils.genomics_utils import opposite_strand, format_chrom
import HTSeq
import numpy as np
import pandas as pd
import seaborn as sns
import six


def handle_args():
    parser = ReferenceArgumentParser(description='Collect stats on read length and biotype')
    parser.add_argument("--intervals", help='.bed/.gtf etc file')
    parser.add_argument("--intervals-name", help="Used in graphs")
    parser.add_argument("--reverse-strand", action='store_true',
                        help="Reverse strand before testing region, useful when you have stranded sequencing and "
                             "the read sequenced is anti-sense")
    parser.add_argument("bam")
    return parser.parse_args()


def get_counts_by_length(bam, regions_array, has_chr, reverse_strand):
    """ bam: bam file path
        regions_array: genomic array of biotypes which is the output from create_biotype_regions_array
        has_chr: Boolean, do reference chromosome names have "chr"? Can be obtained using reference.has_chr
        reverse_strand: switch the strand of the alignment before counting reads
        Returns: pandas dataframe of counts for each biotype for each read length."""
    
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
    
    df = pd.DataFrame(length_counters, columns=sorted(list(length_counters)))
    df = df.fillna(0).astype(int)
    df = df.sort_index().T
    return df


def create_biotype_regions_array(reference, interesting_biotypes=None):
    """ Reference: reference object,
        interesting_biotypes : List of Strings corresponding to biotype keys (everything else is 'other') """

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
            if gene.name and gene.name.startswith("RNY"):
                return "yRNA"
        return gene.biotype

    regions = HTSeq.GenomicArray("auto", stranded=True, typecode='O')
    for transcript in six.itervalues(reference.transcripts):
        # Antisense: Read is in the region of a transcript, but on the opposite strand.
        antisense_iv = transcript.iv.copy()
        antisense_iv.strand = opposite_strand(antisense_iv.strand)
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


def add_intervals_to_regions_array(regions_array, intervals, intervals_name, has_chr):
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
    
    reference = args.reference
    print("Reference is", reference) #To confirm the annotation you're using is what you intended. 
    
    regions_array = create_biotype_regions_array(reference)

    if args.intervals:
        add_intervals_to_regions_array(regions_array, args.intervals, args.intervals_name, reference.has_chr)

    print("Loaded reference - reading BAM")
    df = get_counts_by_length(args.bam, regions_array, reference.has_chr, args.reverse_strand)
    
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
    
    # Add empty columns (biotypes) for those which had zero counts
    df[sorted(set(biotype_colors.keys()).difference(df.columns))] = 0
    
    # Add empty rows (read lengths) for those which had zero counts
    smallest = min(df.index)
    largest = max(df.index)
    all_read_lengths = range(smallest, largest + 1)
    missing_read_lengths = (sorted(set(all_read_lengths).difference(df.index)))
    if missing_read_lengths:
        missing_df = pd.DataFrame(index=missing_read_lengths, dtype=int, columns=df.columns, data=0)
        df = pd.concat([df, missing_df])

    df = df.sort_index()
    df.to_csv(csv_file)
    
    # Graph data
    labels = sorted(biotype_colors.keys())
    colors = []
    for k in labels:
        colors.append(biotype_colors[k])
    
    sns.set_theme(context='paper', style="ticks", font_scale=1.1)

    print("Total read counts:")
    print(df.sum(axis=0))  # A summary of total counts for all read lengths.
    
    fig = Figure(dpi=300, figsize=(4.8, 3.1))
    
    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(111)
     
    # Make stacked bar chart
    bottom = np.zeros(len(df.index), dtype='i')
    for label in df.columns:
        counts = df[label]
        color = biotype_colors[label]
        _ = ax.bar(df.index, counts, label=label, color=color, bottom=bottom, linewidth=0)
        bottom += counts
    
    # Format chart
    ax.set_xlabel("Length (nt)")
    ax.set_ylabel("Read counts")
    
    _, ymax = ax.get_ylim()
    ax.set_ylim(ymin=0, ymax=ymax*1.02)  # Move maximum slightly above highest bar
    ax.set_xlim(xmin=(min(df.index) - 0.7), xmax=(max(df.index) + 0.7))
    
    # Shrink to fit legend
    fig.tight_layout(rect=[0, 0, 0.7, 1])  # left, bottom, right, top
     
    ax.legend(loc='center left', prop={'size': 8.5}, bbox_to_anchor=(1.01, 0.5))
    mk_path_for_file(graph_image)
     
    canvas = FigureCanvasAgg(fig)
    canvas.print_png(graph_image)


if __name__ == '__main__':
    main()
