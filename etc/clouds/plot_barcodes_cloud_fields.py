#!/usr/bin/env python

from homology.util.extract import persistence_intervals_to_array
from homology.util.visualize import plot_barcodes, plot_barcodes_h0_h1
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import os
import numpy as np




CLI = argparse.ArgumentParser(description='Generate and save barcode plots.')
CLI.add_argument(
    'data',
    help='Directory with model data files')
CLI.add_argument(
    'barcodes',
    help='Directory containing barcode output from Ripser')
CLI.add_argument(
    'out',
    help='Directory to store the generated plots in')
CLI.add_argument(
    '--format',
    help='Format for plots to be saved in, e.g. pdf or png',
    default='png')
CLI.add_argument(
    '--sim',
    help='List of simulation names to generate the plots for',
    nargs='*',
    dest='simulations',
    default=['20130401_imicro2', '20130716_imicro2', '20130717_imicro2', '20140717_imicro2'])


args = CLI.parse_args()
data_dir = args.data
barcodes_top_dir = args.barcodes
out_dir = args.out
simulations = args.simulations
plot_format = args.format

for simulation in simulations:
    print 'Processing data from %s' % simulation
    dataset = Dataset('%s/%s/fielddump.ql.001.nc' % (data_dir, simulation))

    barcodes_dir = '%s/%s' % (barcodes_top_dir, simulation)
    filelist = [f for f in os.listdir(barcodes_dir)]
    
    for f in filelist:
        t = int(f.split('_t')[1][:-4])
        print 'timestep %d' % t

        outfile = '%s/%s/clouds_barcode_t%d.%s' % (out_dir, simulation, t, plot_format)
        ql_max = np.amax(dataset.variables['ql'][t,:], 0)
        intervals = persistence_intervals_to_array('%s/%s' % (barcodes_dir, f))
        # Check if return value is None, this indicates no persistence intervals were present
        if intervals is None:
            continue
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25, 25))
        fig.suptitle('Simulation %s, timestep %d' % (simulation, t), fontsize=16)
        ax1.imshow(ql_max, cmap='Blues')
        ax1.invert_yaxis()
        ax1.set_title('QL maximum values')
        # Check if there is actually any H1 homology, otherwise call the simpler function with only the H0 barcodes
        if len(intervals.values()) > 1 and len(intervals.values()[1]) > 0:
            plot_barcodes_h0_h1(ax2, intervals)
        else:
            plot_barcodes(ax2, intervals.values()[0])
        ax2.set_aspect(aspect='auto')
        ax2.set_title('Persistence intervals')
        fig.tight_layout()
        fig.subplots_adjust(top=0.94)
        plt.savefig(outfile, bbox_inches='tight')
        plt.close(fig)

    dataset.close()