#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import subprocess
import os
import argparse

CLI = argparse.ArgumentParser(
    description='Compute persistence intervals for a given set of point clouds.')
CLI.add_argument(
    'topdir',
    help='Top level directory for persistent homology computations.')
CLI.add_argument(
    '--sim',
    help='List of simulation names for which to generate the plots',
    nargs='*',
    dest='simulations',
    default=['20130401_imicro2', '20130716_imicro2', '20130717_imicro2', '20140717_imicro2'])


def main(args):
    top_dir = args.topdir
    simulations = args.simulations

    for simulation in simulations:
        print 'Processing simulation %s ...' % simulation
        # Directory with point cloud files
        points_dir = '%s/%s' % (top_dir, simulation)
        # Where to save the PI files
        out_dir = '%s/barcodes/%s' % (top_dir, simulation)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        # Location of Ripser binary
        ripser_bin = '/usr/local/bin/ripser/ripser'

        filelist = [f for f in os.listdir(
            points_dir) if f.endswith('physicalpoints.csv')]
        for f in filelist:

            fname_in = os.path.join(points_dir, f)
            # print fname_in

            fname_out = '%s/intervals_%s.csv' % (out_dir,
                                                 f.split('_physicalpoints')[0])
            # print fname_out

            with open(fname_out, 'w') as outfile:
                args = ('%s --format point-cloud %s' %
                        (ripser_bin, fname_in)).split()
                p = subprocess.Popen(args, stdout=outfile).communicate()

main(CLI.parse_args())
