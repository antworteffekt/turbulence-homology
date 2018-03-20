#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import AxesGrid
from collections import OrderedDict

def persistence_diagram(file_name, n_steps=None):

    input_file = open(file_name, 'r')
    intervals = input_file.read()
    input_file.close()

    intervals = intervals.split('\n')
    # check if last element is empty, could be EOF
    if len(intervals[-1]) == 0:
        del intervals[-1]
    intervals = [x.split(' ') for x in intervals]
    births = np.array([(int(x[0])) for x in intervals])
    deaths = np.array([(int(x[1])) for x in intervals])
    min_birth = min(births)
    max_death = max(deaths)

    if max_death < 0:
        max_death = max(births) + 1
    if min_birth == max_death:
        max_death = max_death + 1

    # find indices of intervals that die
    finite_intervals = np.where(deaths != -1)

    fig = plt.figure()
    
    plt.plot(births[finite_intervals], deaths[finite_intervals], 'bo')
    # the following is perhaps unnecessary:
    # plt.axis([min_birth, max_death + (max_death - min_birth)/20, min_birth, max_death + (max_death - min_birth)/20])

    # indices of intervals that persist through the entire filtration
    inf_intervals = np.where(deaths == -1)
    inf_vec = (max_death + (max_death - min_birth)/20) * np.ones(len(inf_intervals[0]))
    if len(births[inf_intervals] > 0):
        plt.plot(births[inf_intervals], inf_vec, 'rD')

    # plot the diagonal
    d = range(min_birth, max_death, max(1, (max_death - min_birth) / 20))
    plt.plot(d, d, 'g')
    plt.xlabel('birth')
    plt.ylabel('death')
    # plt.show()

    return fig

def plot_barcodes(ax, intervals, param_dict={}):
    # Requires an array as produced by extract.persistence_intervals_to_array
    # ax : Axes  ---  An axes object to draw to
    ax.set_ylim((-5, intervals.shape[0] + 5))
    ax.set_xlim((0, np.amax(intervals)))
    ax.set_title('Persistence intervals')
    i = 0
    for interval in intervals:
        out = ax.plot((interval[0], interval[1]), (i, i), 'b-', **param_dict)
        i += 1
    return out

def plot_barcodes_h0_h1(ax, intervals, param_dict={}):
    # Requires a dict of arrays as produced by extract.persistence_intervals_to_array
    colors = ['tab:green', 'tab:red', 'tab:purple', 'tab:blue', 'tab:orange']

    n_intervals = sum([x.shape[0] for x in intervals.values()])
    max_value = max([np.amax(x) for x in intervals.values()])

    ax.set_ylim((-5, n_intervals + 1))
    ax.set_xlim((0, max_value))

    idx = 0
    c = 0
    for k, d in intervals.iteritems():
        identifier = k.split('intervals in ')[1]
        for interval in d:
            # print interval
            out = ax.plot((interval[0], interval[1]), (idx, idx), color=colors[c], linestyle='solid',
                           label=identifier)
            idx += 1
        c += 1
    # The following lines remove all duplicate values in handles and labels before passing them to legend()
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())
    return out


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap
