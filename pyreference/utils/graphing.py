'''
Created on 22Jan.,2018

@author: dlawrence
'''

import abc
import logging
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.cm import ScalarMappable
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import os

import numpy as np


def axis_stacked_bar_graph(ax, largs, arrays, labels, colors):
    bottom = np.zeros(len(largs), dtype='i')
    for array, label, color in zip(arrays, labels, colors):
        lines = ax.bar(largs, array, label=label, color=color, bottom=bottom, linewidth=0)
        bottom += array
    
    return lines


def write_stacked_bar_graph(graph_image, largs, arrays, labels, colors, **kwargs):
    '''largs = x-values
    arrays = list of lists y-values
    labels = list of labels, same length as arrays
    colors = list of colors, same length as arrays

    kwargs:
        extra_func = method called with args=(ax, fig, lines)
        legend_side_ratio : eg 0.1 - Shrink figure by this
        legend_kwargs : Passed to ax.legend
    '''
    logging.info("write_stacked_bar_graph: %s", graph_image)

    title = kwargs.get("title")
    extra_func = kwargs.get("extra_func")
    subtitle = kwargs.get("subtitle")
    x_label = kwargs.get("x_label")
    y_label = kwargs.get("y_label")
    legend_side_ratio = kwargs.get("legend_side_ratio")

    # Old kwargs, which we passed to legend
    legend_kwargs = {"loc" : kwargs.get("loc"),
                     "title" : kwargs.get("legend_title"),
                     "prop" : kwargs.get("legend_props", {})}
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

    ax.set_xlim(xmin=largs[0], xmax=largs[-1]+1)
    
    if legend_side_ratio is not None:
        print("legend_side_ratio: %f" % legend_side_ratio)
        assert 0 < legend_side_ratio < 1, "legend_side_ratio: 0 < %f < 1" % legend_side_ratio 
        
        # Shrink current axis's height by legend_below_ratio on the bottom
        #ax = pyplot.axes()
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
    ''' same as write_stacked_bar_graph but writes out 1 plot per entry in array '''

    (name_base, extension) = os.path.splitext(graph_image)
    
    title = kwargs.get("title", "")
    y_label = kwargs.get("y_label", None)
    for array, label, color in zip(arrays, labels, colors):
        kwargs["title"] = "%s (%s)" % (title, label)
        kwargs["color"] = color

        bargraph = BarGraph(largs, array)
        bargraph.y_label = y_label
        individual_graph_image = "%s.%s%s" % (name_base, label, extension) # extension still has dot
        bargraph.save(individual_graph_image)
    
    
class GraphBase(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, **kwargs):
        '''
            legend: [('label', 'color'), ('label2', 'color2')]
        '''
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
        ''' a hook method if you want to do something about the figure '''
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
            plot(ax) # Implementation
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
    
    
    
class BarGraph(GraphBase):
    def __init__(self, x, y, **kwargs):
        super(BarGraph, self).__init__(**kwargs)
        self.x = x
        self.y = y
        self.labels = kwargs.get("labels")
        self.color = kwargs.get("color")
        
    def plot(self, ax):
        width = 0.5

        ax.set_xlim(xmin=self.x[0], xmax=self.x[-1]+1)
        ax.bar(self.x, self.y, width=width, color=self.color)

        if self.labels:
            ax.set_xticks(self.x + width/2)
            ax.set_xticklabels(self.labels)