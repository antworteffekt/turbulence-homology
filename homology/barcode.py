# -*- coding: utf-8 -*-
"""
Created on Tue May 29 15:45:23 2018

@author: licon
"""

from __future__ import print_function
from homology.interval import Interval
import re
from collections import OrderedDict
import numpy as np
import subprocess
import StringIO

class Barcode:

    def __init__(self, value_range=None, max_dim=0):
        """
        Generates empty barcode. Will be filled with data separately.
        """
        self.value_range = value_range
        self.max_dim = max_dim
        self.intervals = {}
        self.additive_persistence = {}
        self.multiplicative_persistence = {}

    def __repr__(self):

        barcode_str = 'Persistence barcode:\n\tmaximum dimension:{}\n\tValue range:{}\n'.format(
            self.max_dim,
            self.value_range
        )
        intervals_str = ''
        for k, v in self.intervals.items():
            intervals_str += 'Dimension {} : {} intervals\n'.format(k, len(v))
        return(barcode_str + intervals_str)


    def __str__(self):
        string_rep = ''
        value_range = str(self.value_range)
        value_range = value_range.replace('(', '[')
        value_range = value_range.replace(')', ']')
        string_rep += "value range: %s" % value_range
        i = 0
        for dim in self.intervals.values():
            string_rep += "\npersistence intervals in dim %d:\n" % i
            i += 1
            # Remove first and last characters
            dim = str(dim)[1:-1]
            dim = dim.replace('(', '[')
            dim = dim.replace(', [', '\n[')
            string_rep += dim
        return string_rep

    def __getitem__(self, dim):
        if dim not in self.intervals.keys():
            raise IndexError(
                'Intervals for dimension {} not found in barcode.'.format(dim))
        return self.intervals[dim]

    def compute(self, data, dim=1, data_format='point-cloud'):
        # Compute and parse in one call
        intervals = self.persistence_intervals(data, '/home/licon/src/ripser/ripser', dim, data_format)
        self.parse_ripser_output(intervals)

    def persistence_intervals(self, X, ripser_bin, dim, data_format='point-cloud'):
    
        x_str = StringIO.StringIO()
        np.savetxt(x_str, X, delimiter=',')
        # Validate data format
        if not data_format in ['point-cloud', 'distance']:
            raise ValueError("Data format not understood: %s" % data_format)
        args = [ripser_bin, '--format', data_format, '--dim', str(dim)]
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE, 
                             stderr=subprocess.STDOUT)
        ripser_out = p.communicate(input=x_str.getvalue())[0]
        return ripser_out

    def parse_ripser_output(self, barcode_str):

        barcode_str = barcode_str.split('\n')
        if barcode_str[2].startswith('value'):
            # Assume the value range is given in position 2
            value_range = barcode_str[2].split(': ')[1]
        elif barcode_str[1].startswith('value'):
            value_range = barcode_str[1].split(': ')[1]
        # Remove brackets
        value_range = value_range[1:-1]
        value_range = (float(value_range.split(',')[0]), float(
            value_range.split(',')[1]))
        self.value_range = value_range

        dim_strings = [s for s in barcode_str if s.startswith(
            'persistence intervals in dim')]

        dimensions = [int(re.sub('[^0-9]', '', d)) for d in dim_strings]
        self.max_dim = max(dimensions)

        self.intervals = {x: [] for x in range(self.max_dim + 1)}
        dim_idx = [barcode_str.index(s) for s in dim_strings]
        dim_idx.append(-1)

        for i in range(1, len(dim_idx)):

            intervals = barcode_str[dim_idx[i - 1]:dim_idx[i]]
            identifier = intervals.pop(0).replace(':', '')[-1]
            identifier = int(identifier)
            inf_intervals = [s for s in intervals if s.endswith(', )')]
            intervals = [s for s in intervals if not s.endswith(', )')]

            # Remove brackets
            intervals = map(lambda x: x.replace('[', ''), intervals)
            intervals = map(lambda x: x[:-1], intervals)
            intervals = [(float(s.split(',')[0]), float(s.split(',')[1]))
                         for s in intervals]
            # Add infinite intervals, if present:
            if inf_intervals:
                inf_intervals = map(
                    lambda x: x.replace('[', ''), inf_intervals)
                inf_intervals = map(lambda x: x[:-1], inf_intervals)
                inf_intervals = [(float(s.split(',')[0]), )
                                 for s in inf_intervals]
                intervals = inf_intervals + intervals
            self.intervals[identifier] = intervals

    def plot(self, ax, dimensions=None, plot_infinite_intervals=False):
        # Requires a dict of lists as produced by extract.parse_output
        # colors = ['tab:red', 'tab:blue',
        #           'tab:purple', 'tab:green', 'tab:orange']
        colors = {0: 'forestgreen',
                  1: 'red',
                  2: 'tab:purple',
                  3: 'tab:green',
                  4: 'tab:orange'}

        if plot_infinite_intervals:
            ax.set_xlim([self.value_range[0], self.value_range[1]])
        else:
            ends = [x for sublist in self.intervals.values() for x in sublist]
            ends = [x[1] for x in ends if len(x) > 1]
            ax.set_xlim([self.value_range[0], max(ends) * 1.05])

        idx = 0
        c = 0
        if dimensions is None:
            # Plot intervals in all dimensions
            n_intervals = sum([len(v) for k, v in self.intervals.items()])
            ax.set_ylim([-5, n_intervals + 1])
            for k, d in self.intervals.iteritems():
                identifier = 'dim %s' % k
                for interval in d:
                    if len(interval) == 1 and plot_infinite_intervals:
                        interval = (interval[0], self.value_range[1])
                    elif len(interval) == 1:
                        continue
                    out = ax.plot((interval[0], interval[1]), (idx, idx), color=colors[k], linestyle='solid',
                                  lw=0.99, alpha=0.9, label=identifier)
                    idx += 1
                c += 1
        else:
            # Plot only the specified dimensions
            subintervals = {k: self.intervals[k] for k in dimensions if k in self.intervals}
            n_intervals = sum([len(v) for k, v in subintervals.items()])
            ax.set_ylim([-5, n_intervals + 1])
            for k, d in subintervals.iteritems():
                identifier = 'dim %s' % k
                for interval in d:
                    if len(interval) == 1 and plot_infinite_intervals:
                        interval = (interval[0], self.value_range[1])
                    elif len(interval) == 1:
                        continue
                    out = ax.plot((interval[0], interval[1]), (idx, idx), color=colors[k], linestyle='solid',
                                  lw=0.99, alpha=0.9, label=identifier)
                    idx += 1
                c += 1
        # The following lines remove all duplicate values in handles and labels
        # before passing them to legend()
        handles, labels = ax.get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
        ax.get_yaxis().set_visible(False)
        return out

    def to_string(self):
        # if barcode is None:
        #     return 'None'
        string_rep = ''
        value_range = str(self.value_range)
        value_range = value_range.replace('(', '[')
        value_range = value_range.replace(')', ']')
        string_rep += "value range: %s" % value_range
        i = 0
        for dim in self.intervals:
            string_rep += "\npersistence intervals in dim %d:\n" % i
            i += 1
            # Remove first and last characters
            dim = str(dim)[1:-1]
            dim = dim.replace('(', '[')
            dim = dim.replace(', [', '\n[')
            string_rep += dim
        return string_rep

    def S_a(self, xs, dim):
        """
        Evaluate the feature-counting invariant for additive persistence.
        """
        # If persistence values have not been precomputed, compute them now. Sort the resulting array.
        try:
            if not bool(self.additive_persistence):
                self.additive_persistence = {
                    k : np.sort(np.array([x[1] - x[0] for x in v if len(x) > 1])) for k, v in self.intervals.iteritems()
                }
        except AttributeError:
            self.additive_persistence = {
                    k : np.sort(np.array([x[1] - x[0] for x in v if len(x) > 1])) for k, v in self.intervals.iteritems()
                }
        ys = np.sum(self.additive_persistence[dim] > xs[:, None], axis=1)
        return ys

    def S_p(self, xs, dim):
        """
        Evaluate the feature-counting invariant for multiplicative persistence.
        """
        # If persistence values have not been precomputed, compute them now. Sort the resulting array.
        try:
            if not bool(self.multiplicative_persistence):
                self.multiplicative_persistence = {
                    k : [x[1] / x[0] for x in v if x[0] > 1e-10] for k, v in self.intervals.iteritems() if k > 0
                }
        except AttributeError:
            self.multiplicative_persistence = {
                    k : [x[1] / x[0] for x in v if x[0] > 1e-10] for k, v in self.intervals.iteritems() if k > 0
                }

        ys = np.sum(self.multiplicative_persistence[dim] > xs[:, None], axis=1)
        return ys

if __name__ == '__main__':

    import matplotlib.pyplot as plt
    interval = Interval(start=10, end=15)
#    print(interval)

    barcode = Barcode(value_range=(10, 100), max_dim=5)

    infile = ('/home/licon/uni-koeln/persistence/pointcloud_data_DALES/first_barcode_set_21032018/barcodes/20130401_imicro2/PI_clouds_2dsample_t45.csv')

    with open(infile, 'r') as f:
        barcode_str = f.read()

    barcode.parse_ripser_output(barcode_str)
    print(barcode)
    fig, ax = plt.subplots()
    barcode.plot(ax)
    
