#!/usr/bin/env python

from __future__ import print_function
import numpy as np
from netCDF4 import Dataset
import argparse
import os
from homology.sampler import Sampler
import pickle
import config
from itertools import product
from functools import partial
import multiprocessing as mp
from collections import defaultdict

def process_arguments():

    parser = argparse.ArgumentParser(description="Compute size of connected components given threshold value.")

    parser.add_argument(
        'env',
        help='Which model environment to use (one of LES, DNS, DALES, TEST)')

    parser.add_argument(
        '-v', '--varname',
        help='Variable name. {x,y,z}_wind',
        required=True)

    parser.add_argument(
        '-s',
        help='Simulation name. For the LES case, SPx. For DALES, it is of the form 20140717_imicro2.',
        required=False,
        default='')

    parser.add_argument(
        '--multiprocess',
        help='Whether to use the multiprocessing version.',
        action='store_true')

    parser.add_argument(
        '--nproc',
        help='Number of processors to use, if applicable.',
        default=2,
        type=int)

    return parser.parse_args()


def structure_sizes_multiprocess(input_fname, varname, rescale, axis, periodicity, params):
# def structure_sizes_multiprocess(input_fname, varname, rescale, axis, params):    
    """
    Given a dataset, rescale (if necessary), split into connected components,
    and aggregate their sizes. Return as a list of integers. 
    IN
        input_fname  : path to netCDF file containing the desired data
        varname      : variable name
        axis         : position of the spatial axis to iterate over
        periodicity  :
        params       : list of [time, height, threshold]
        **kwargs     : other parameters; rescale? If true, then it also contains 
                       offset and scale values.
    OUT
        list of integers representing the sizes of the connected components.

    """
    # Unpack parameters
    t, v, threshold = params
    # print ('Timestep %d, Height %d' % (t, v))
    # Open dataset connection
    dataset = Dataset(input_fname, 'r')

    idx = [t, slice(None, None, None), slice(None, None, None), slice(None, None, None)]
    idx[axis] = v
    # print(idx)
    data = dataset.variables[varname][idx]
    # print (data.shape)
    if rescale:
        data = dataset.variables[varname].var_add_offset + \
               dataset.variables[varname].var_scale_factor * data

    # Close dataset connection
    dataset.close()

    if threshold == 'mean':
        sampler = Sampler(data > np.mean(data))
        sampler.connected_components(periodic=periodicity)
    else:
        sampler = Sampler(data > threshold)
        sampler.connected_components(periodic=periodicity)

    sizes = [x for x in sampler.uf.size if x != 0]
    # print (sizes)
    return {'time' : t,
            'var' : v,
            'threshold' : threshold,
            'sizes' : sizes}


if __name__ == '__main__':
    """
    Load file, separate components, calculate sizes.
    """

    # TODO: add boolean flag to generate test data only (2 timesteps, two plane cross sections)

    args = process_arguments()
    # print (args)

    # Select environment
    if args.env == 'LES':
        ci = config.LES
    elif args.env == 'DNS':
        ci = config.DNS
    elif args.env == 'DALES':
        ci = config.DALES
    elif args.env == 'TEST':
        ci = config.TEST
    else:
        raise ValueError("Invalid model environment selected.")
    
    # Variable name
    var_key = args.varname
    var_name = ci.VARIABLES[var_key]
    # Choose directions to do the array slicing
    if var_key == 'x_wind':
        dirs = ci.DIRECTIONS['x']
    elif var_key == 'y_wind':
        dirs = ci.DIRECTIONS['y']
    elif var_key == 'z_wind':
        dirs = ci.DIRECTIONS['z']

    simulation = args.s
    # Output directory. Need simulation appended
    out_dir = '%s/%s' % (ci.OUT, simulation)
    # print (out_dir)
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # Check if filename contains variable name
    if 'VARIABLE' in ci.FILENAME:
        ci.FILENAME = ci.FILENAME.replace('VARIABLE', ci.VARIABLES[var_key])
    
    # Check if size needs to be included as well - this should only be the case for DNS dataset
    if 'SIZE' in ci.FILENAME:
        ci.FILENAME = ci.FILENAME.replace('SIZE', ci.SIZES[simulation])
    # Append simulation name as well.
    in_fname = '%s/%s/%s' % (ci.FILES, simulation, ci.FILENAME)
    
    thresholds = ci.THRESHOLDS

    # print(in_fname)
    # print(out_dir)

    # Open file connection to read parameters
    d = Dataset(in_fname, 'r')
    # Time range
    t_range = range(d.dimensions[ci.DIMENSIONS['t']].size)[50:51]

    # Variable range. Given by the first character in the case of velocity vector component (z_wind, etc.)
    var_direction = var_key[0]
    dim_name = ci.DIMENSIONS[var_direction]
    try:
        axis = d[var_name].dimensions.index(dim_name)
    except ValueError:
        print ("Dimension name %s not found in dataset %s for variable %s. Trying alternate names ..." % 
            (dim_name, in_fname, var_name))

        if '%s_stag' % dim_name in d[var_name].dimensions:
            dim_name = '%s_stag' % dim_name
            axis = d[var_name].dimensions.index(dim_name)
        elif '%sm' % var_direction in d[var_name].dimensions:
            dim_name = '%sm' % var_direction
            axis = d[var_name].dimensions.index(dim_name)
        else:
            raise
        
    var_range = range(d.dimensions[dim_name].size)
    d.close()
    # exit()

    if args.multiprocess:

        # print ("Multiprocess")

        func = partial(structure_sizes_multiprocess, in_fname, var_name, ci.RESCALE, axis, ci.PERIODIC)
        params = product(t_range, var_range, thresholds.values())
        
        pool = mp.Pool(processes=args.nproc)
        # for p in params:
        #     print(p)
        res = pool.map(func, params)
        pool.close()

        # Reorganize data for writing
        sizes = defaultdict(dict)
        for t, v in product(t_range, var_range):
            sizes[t][v] = [x['sizes'] for x in res if x['time'] == t and x['var'] == v][0]

        # TODO :change the way threshold values are handled
        k = thresholds.keys()[0]
        out_fname = '%s/%s_structure_sizes_%s.pickle' % (out_dir, var_key, k)
        with open(out_fname, 'w') as f:
            pickle.dump(sizes, f)

    else:

        d = Dataset(in_fname, 'r')
        # print (in_fname)
        for k, threshold in thresholds.items():

            out_fname = '%s/%s_structure_sizes_%s.pickle' % (out_dir, var_key, k)
            sizes = {}

            for t in t_range:
                sizes[t] = {}
                for v in var_range:
                    # Indices for slicing
                    idx = [t, slice(None, None, None), slice(None, None, None), slice(None, None, None)]
                    idx[axis] = v
                    # print(idx)
                    # separate components in two-dimensional field
                    X = d.variables[var_name][idx]
                    if ci.RESCALE:
                        # print("rescaling")
                        X = d[var_name].var_add_offset + d[var_name].var_scale_factor * X
                    sampler = Sampler(X > threshold)
                    # print (X.shape)
                    # print (np.mean(X))
                    # print("periodicity: %s" % ci.PERIODIC)
                    # print("threshold: %s" % threshold)
                    sampler.connected_components(periodic=ci.PERIODIC)


                    # print([x for x in sampler.uf.size if x != 0])

                    sizes[t][v] = [x for x in sampler.uf.size if x != 0]

            print ("Done.")
        d.close()

        # # Check if output file exists
        # if os.path.isfile(out_fname):
        #     # Open as read + write
        #     with open(out_fname, 'r') as f:
        #         sizes_old = pickle.load(f)
            
        #     sizes.update(sizes_old)

        #     with open(out_fname, 'w') as f:
        #         pickle.dump(sizes, f)

        # else:
        #     # Open as write
        with open(out_fname, 'w') as f:
            pickle.dump(sizes, f)

    # for k, d in sizes.items():
    #     print ("**********\n   Timestep %d :"%k)
    #     for u, v in d.items():
    #         print ("level %d : %d components\n%s" % (u, len(v), v))



