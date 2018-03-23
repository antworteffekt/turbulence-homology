#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import argparse
import os
from homology.util.sample import *

def process_arguments():

    parser = argparse.ArgumentParser(
        description="""
        Generate point samples to be used as point cloud input for persistent homology.
        Will use the ~/var/tmp/ directory to store temporary output files; it will be created if it 
        does not exist. 

        A given number n of samples per simulation timestep will be created.
        """)

    parser.add_argument(
        '-i', '--input-file',
        help='Path to netCDF file containing the data to be sampled.',
        required=True)
    parser.add_argument(
        '-o', '--output-dir',
        help='Path to top level directory where samples will be stored.',
        required=True)
    parser.add_argument(
        '--varname',
        help='Variable name to use for sampling.',
        required=True)
    parser.add_argument(
        '--project',
        help='Whether to use the column-wise projection of the selected variable.',
        action='store_true')
    parser.add_argument(
        '--sample-ratio',
        help='Proportion of domain size to sample.',
        type=float,
        default=0.05)
    parser.add_argument(
        '-n',
        help='Number of samples to generate.',
        default=10,
        type=int)
    
    try:
        return parser.parse_args()
    except IOError, msg:
        parser.error(str(msg))


def main():
    
    args = process_arguments()
    # print args
    # Seed the generator with a random prime number, say 57
    np.random.seed(57)
    # Check for ~/var/tmp
    tmp = os.path.expanduser('~/var/tmp')
    if not os.path.exists(tmp):
        os.makedirs(tmp)

    # Check for top level output directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Open file connection
    dataset = Dataset(args.input_file, 'r')
    timerange = range(dataset[args.varname].shape[0])
    
    # Iterate over simulation timesteps
    for t in timerange:
        print '********* Timestep: %d *********' % t
        if args.project:
            # Get two-dimensional projected field
            domain = project_2d_cloud_field(dataset, t, args.varname)

            # Get the required number of samples and write to disk
            for i in range(args.n):
                sample = sample_connected_components(domain, tmp, args.sample_ratio)
                if sample is None:
                    continue
                else:
                    out_fname = '%s/t%d_sample%d' % (args.output_dir, t, i+1)
                    save_sample(dataset, sample, out_fname)
        else:
            # Assume the samples will be taken in horizontal planes
            assert dataset.variables[args.varname].dimensions[1] == 'zm'
            z_range = range(dataset.variables[args.varname].shape[1])
            for z in z_range:
                # Create separate directories for each horizontal segment
                current_out_dir = '%s/z%d' % (args.output_dir, z)
                if not os.path.exists(current_out_dir):
                    os.makedirs(current_out_dir)
                domain = (dataset.variables[args.varname][t,z,:] > 0)
                # Get the required number of samples and write to disk
                for i in range(args.n):
                    sample = sample_connected_components(domain, tmp, args.sample_ratio)
                    if sample is None:
                        continue
                    else:
                        out_fname = '%s/t%d_sample%d' % (current_out_dir, t, i+1)
                        save_sample(dataset, sample, out_fname)

    dataset.close()

if __name__ == '__main__':
    main()