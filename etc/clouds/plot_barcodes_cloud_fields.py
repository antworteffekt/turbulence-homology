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
    'topdir',
    help='Top level directory for persistent homology computations')
CLI.add_argument(
    '--format',
    help='Format for plots to be saved in, e.g. pdf or png',
    default='png')
CLI.add_argument(
    '--sim',
    help='List of simulation names for which to generate the plots',
    nargs='*',
    dest='simulations',
    default=['20130401_imicro2', '20130716_imicro2', '20130717_imicro2', '20140717_imicro2'])
CLI.add_argument(
    '--binary',
    help='Whether to plot the cloud fields as binary domains',
    action='store_true')

def main(args):
    data_dir = args.data
    top_dir = args.topdir
    simulations = args.simulations
    plot_format = args.format
    plot_binary_field = args.binary

    out_dir = '%s/plots' % top_dir

    for simulation in simulations:
        print 'Processing data from %s' % simulation
        dataset = Dataset('%s/%s/fielddump.ql.001.nc' % (data_dir, simulation))

        barcodes_dir = '%s/barcodes/%s' % (top_dir, simulation)
        filelist = [f for f in os.listdir(barcodes_dir)]
        
        for f in filelist:
            t = int(f.split('_t')[1][:-4])
            print 'timestep %d' % t
            # Compute maximum values
            ql_max = np.amax(dataset.variables['ql'][t,:], 0)
            # Load cloud sample points into an array
            sample_fname = '%s/samples/%s/gridpoints_2dsample_t%d.csv' % (top_dir, simulation, t)
            sampled_points = np.loadtxt(sample_fname, delimiter=',')
            # Convert saved persistence intervals to numpy array
            intervals = persistence_intervals_to_array('%s/%s' % (barcodes_dir, f))
            # Check if return value is None, this indicates no persistence intervals were present
            if intervals is None:
                continue

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25, 25))
            fig.suptitle('Simulation %s, timestep %d' % (simulation, t), fontsize=16)
            if plot_binary_field:
                ax1.imshow(ql_max > 0, cmap='Greys', interpolation='none')
                ax1.scatter(sampled_points[:,1], sampled_points[:,0],
                            facecolors='none',
                            edgecolors='aqua',
                            alpha=1,
                            marker='o')

            else:
                ax1.imshow(ql_max, cmap='Blues')
                ax1.scatter(sampled_points[:,1], sampled_points[:,0],
                            facecolors='none',
                            edgecolors='crimson',
                            alpha=0.6,
                            marker='o')

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
            # Save plot to disk
            outfile = '%s/%s/clouds_barcode_t%d.%s' % (out_dir, simulation, t, plot_format)
            plt.savefig(outfile, bbox_inches='tight')
            plt.close(fig)

        dataset.close()

main(CLI.parse_args())