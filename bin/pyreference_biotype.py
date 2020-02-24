from __future__ import print_function, absolute_import

import HTSeq
import abc
from argparse import ArgumentParser
from collections import Counter
import csv
import logging
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.cm import ScalarMappable
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import os
from pyreference import Reference
from pyreference.utils import iv_iterators
from pyreference.utils.file_utils import name_from_file_name, file_or_file_name
from pyreference.utils.genomics_utils import opposite_strand, format_chrom
import six
import sys

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


def create_biotype_regions_array(reference, interesting_biotypes=None):
    """ genes_by_biotype : dict of {"biotype" : genes[]}
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
            if gene.name.startswith("RNY"):
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

    regions_array = create_biotype_regions_array(reference)

    if args.intervals:
        add_intervals_to_regions_array(regions_array, args.intervals, args.intervals_name, reference.has_chr,
                                       args.reverse_strand)

    print("Loaded reference - reading BAM")
    length_counters = get_length_counts(args.bam, regions_array, reference.has_chr, args.reverse_strand)

    start = six.MAXSIZE
    end = 0
    for (i, c) in enumerate(length_counters):
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

    graph_types = {"read_counts": lambda x: x}

    sample_name = name_from_file_name(args.bam)

    for (y_label, op) in six.iteritems(graph_types):
        graph_image = "%s.%s.regions.png" % (sample_name, y_label)
        data = [op(i) for i in arrays]

        ncol = max(len(labels), 5)
        legend_kwargs = {'loc': 'center left',
                         'prop': {'size': 8},
                         'bbox_to_anchor': (1.01, 0.5)}
        write_stacked_bar_graph(graph_image, largs, data, labels, colors, y_label=y_label, legend_side_ratio=0.15,
                                legend_kwargs=legend_kwargs)
        write_individual_bar_graphs(graph_image, largs, data, labels, colors, y_label=y_label)

    # Write region stats CSV
    csv_file = "%s.biotype_stats.csv" % sample_name
    write_csv_dict(csv_file, labels, [total_counts])


"""
Created on 22Jan.,2018

@author: dlawrence
"""


def axis_stacked_bar_graph(ax, largs, arrays, labels, colors):
    bottom = np.zeros(len(largs), dtype='i')
    lines = None
    for array, label, color in zip(arrays, labels, colors):
        lines = ax.bar(largs, array, label=label, color=color, bottom=bottom, linewidth=0)
        bottom += array

    return lines


def write_stacked_bar_graph(graph_image, largs, arrays, labels, colors, **kwargs):
    """largs = x-values
    arrays = list of lists y-values
    labels = list of labels, same length as arrays
    colors = list of colors, same length as arrays

    kwargs:
        extra_func = method called with args=(ax, fig, lines)
        legend_side_ratio : eg 0.1 - Shrink figure by this
        legend_kwargs : Passed to ax.legend
    """
    logging.info("write_stacked_bar_graph: %s", graph_image)

    title = kwargs.get("title")
    extra_func = kwargs.get("extra_func")
    subtitle = kwargs.get("subtitle")
    x_label = kwargs.get("x_label")
    y_label = kwargs.get("y_label")
    legend_side_ratio = kwargs.get("legend_side_ratio")

    # Old kwargs, which we passed to legend
    legend_kwargs = {"loc": kwargs.get("loc"),
                     "title": kwargs.get("legend_title"),
                     "prop": kwargs.get("legend_props", {})}
    # Overwrite with legend_kwargs
    legend_kwargs.update(kwargs.get("legend_kwargs", {}))

    fig = Figure(dpi=300)
    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(111)
    lines = axis_stacked_bar_graph(ax, largs, arrays, labels, colors)

    if title:
        fig.suptitle(title, fontsize=18)
    if subtitle:
        ax.set_title(subtitle, fontsize=12)
    if extra_func:
        extra_func(ax, fig, lines)
    if x_label:
        ax.set_xlabel(x_label)
    if y_label:
        ax.set_ylabel(y_label)

    ax.set_xlim(xmin=largs[0], xmax=largs[-1] + 1)

    if legend_side_ratio is not None:
        print("legend_side_ratio: %f" % legend_side_ratio)
        assert 0 < legend_side_ratio < 1, "legend_side_ratio: 0 < %f < 1" % legend_side_ratio

        # Shrink current axis's height by legend_below_ratio on the bottom
        # ax = pyplot.axes()
        box = ax.get_position()
        print("Box w/h = (%f,%f)" % (box.width, box.height))

        new_width = box.height * (1.0 - legend_side_ratio)
        new_position = [box.x0, box.y0,
                        new_width, box.height]

        ax.set_position(new_position)

        box = ax.get_position()
        print("NEW Box w/h = (%f,%f)" % (box.width, box.height))

    print("**legend_kwargs = %s" % legend_kwargs)
    ax.legend(**legend_kwargs)

    canvas = FigureCanvasAgg(fig)
    canvas.print_png(graph_image)


def write_individual_bar_graphs(graph_image, largs, arrays, labels, colors, **kwargs):
    """ same as write_stacked_bar_graph but writes out 1 plot per entry in array """

    (name_base, extension) = os.path.splitext(graph_image)

    title = kwargs.get("title", "")
    y_label = kwargs.get("y_label", None)
    for array, label, color in zip(arrays, labels, colors):
        kwargs["title"] = "%s (%s)" % (title, label)
        kwargs["color"] = color

        bargraph = BarGraph(largs, array)
        bargraph.y_label = y_label
        individual_graph_image = "%s.%s%s" % (name_base, label, extension)  # extension still has dot
        bargraph.save(individual_graph_image)


class GraphBase(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, **kwargs):
        """
            legend: [('label', 'color'), ('label2', 'color2')]
        """
        self.title = kwargs.get("title")
        self.x_label = kwargs.get("x_label")
        self.y_label = kwargs.get("y_label")
        self.legend = kwargs.get('legend')
        self.plot_methods = [self.plot]

    def decorations(self, ax):
        if self.title:
            ax.set_title(self.title)
        if self.x_label:
            ax.set_xlabel(self.x_label)
        if self.y_label:
            ax.set_ylabel(self.y_label)

    @abc.abstractmethod
    def plot(self, ax):
        return

    def post_plot(self, ax):
        if self.legend:
            patches = []
            labels = []
            for (name, color) in self.legend:
                labels.append(name)
                patches.append(Rectangle((0, 0), 1, 1, fc=color))

            ax.legend(patches, labels, loc='upper left')

    def figure(self, figure):
        """ a hook method if you want to do something about the figure """
        pass

    def get_figure_and_axis(self, dpi):
        figure = Figure(dpi=dpi)
        figure.patch.set_facecolor('white')
        ax = figure.add_subplot(1, 1, 1)
        return figure, ax

    def save(self, filename_or_obj, dpi=None, file_type='png'):
        figure, ax = self.get_figure_and_axis(dpi)

        self.decorations(ax)
        for plot in self.plot_methods:
            plot(ax)  # Implementation
        self.post_plot(ax)
        self.figure(figure)

        canvas = FigureCanvasAgg(figure)
        if file_type == 'png':
            canvas.print_png(filename_or_obj)
        elif file_type == 'pdf':
            canvas.print_pdf(filename_or_obj)

    def draw_origin(self, ax):
        ax.axvline(c='black', lw=0.5)
        ax.axhline(c='black', lw=0.5)

    def colorbar_from_cmap_array(self, figure, cmap, array):
        mappable = ScalarMappable(cmap=cmap)
        mappable.set_array(array)
        figure.colorbar(mappable)


def write_csv_dict(csv_file, headers, rows, extrasaction=None, dialect=None):
    """
    default dialect = 'excel', other dialect option: 'excel-tab'
    These are the same optional arguments as csv.DictWriter
    headers=keys for dicts
    rows=list of dicts
    """

    if extrasaction is None:
        extrasaction = "raise"
    if dialect is None:
        dialect = 'excel'

    f = file_or_file_name(csv_file, "wb")

    writer = csv.DictWriter(f, headers, extrasaction=extrasaction, dialect=dialect)
    writer.writerow(dict(zip(headers, headers)))
    writer.writerows(rows)


class BarGraph(GraphBase):
    def __init__(self, x, y, **kwargs):
        super(BarGraph, self).__init__(**kwargs)
        self.x = x
        self.y = y
        self.labels = kwargs.get("labels")
        self.color = kwargs.get("color")

    def plot(self, ax):
        width = 0.5

        ax.set_xlim(xmin=self.x[0], xmax=self.x[-1] + 1)
        ax.bar(self.x, self.y, width=width, color=self.color)

        if self.labels:
            ax.set_xticks(self.x + width / 2)
            ax.set_xticklabels(self.labels)


if __name__ == '__main__':
    main()
