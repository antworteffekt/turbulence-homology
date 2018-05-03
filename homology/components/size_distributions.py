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


def select_data(in_fname, multiprocess, varname, rescale, axis, t, v):

    idx = [t, slice(None, None, None), slice(None, None, None), slice(None, None, None)]
    idx[axis] = v

    # Open dataset connection if a string is given, i.e. if multiprocessing flag is signalled
    if multiprocess:
        dataset = Dataset(in_fname, 'r')
        data = dataset.variables[varname][idx]
        if rescale:
            data = dataset.variables[varname].var_add_offset + \
                   dataset.variables[varname].var_scale_factor * data
        # Close dataset connection
        dataset.close()

    else:
        # The other possibility is that input is already a file connection
        dataset = in_fname
        data = dataset.variables[varname][idx]
        if rescale:
            data = dataset.variables[varname].var_add_offset + \
                   dataset.variables[varname].var_scale_factor * data

    return data


def structure_sizes(input_fname, multiprocess, varname, rescale, axis, periodicity, params):
    """
    Given a dataset, rescale (if necessary), split into connected components,
    and aggregate their sizes. Return as a list of integers. 
    IN
        input_fname  : path to netCDF file containing the desired data, or file connection.
        multiprocess : boolean: if True, multiprocessing is being used and the input will be given as string
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
    
    # idx = [t, slice(None, None, None), slice(None, None, None), slice(None, None, None)]
    # idx[axis] = v

    # # Open dataset connection if a string is given, i.e. if multiprocessing flag is signalled
    # if multiprocess:
    #     dataset = Dataset(input_fname, 'r')
    #     data = dataset.variables[varname][idx]
    #     if rescale:
    #         data = dataset.variables[varname].var_add_offset + \
    #                dataset.variables[varname].var_scale_factor * data
    #     # Close dataset connection
    #     dataset.close()

    # else:
    #     # The other possibility is that input is already a file connection
    #     dataset = input_fname
    #     data = dataset.variables[varname][idx]
    #     if rescale:
    #         data = dataset.variables[varname].var_add_offset + \
    #                dataset.variables[varname].var_scale_factor * data
    X = select_data(input_fname, multiprocess, varname, rescale, axis, t, v)

    if threshold == 'mean':
        sampler = Sampler(X > np.mean(data))
        sampler.connected_components(periodic=periodicity)
    else:
        sampler = Sampler(X > threshold)
        sampler.connected_components(periodic=periodicity)

    sizes = [x for x in sampler.uf.size if x != 0]
    return {'time' : t,
            'var' : v,
            'threshold' : threshold,
            'sizes' : sizes}


if __name__ == '__main__':
    """
    Load file, separate components, calculate sizes.
    """

    # TODO: add boolean flag to generate test data only (2 timesteps, two plane cross sections)
    testdata = True

    args = process_arguments()

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

    # Open file connection to read parameters
    d = Dataset(in_fname, 'r')
    # Time range
    if testdata:
        t_range = range(d.dimensions[ci.DIMENSIONS['t']].size)[10:12]
        var_range = range(d.dimensions[dim_name].size)[5:7]
    else:
        t_range = range(d.dimensions[ci.DIMENSIONS['t']].size)
        var_range = range(d.dimensions[dim_name].size)

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

    d.close()

    # TODO :change the way threshold values are handled
    k = thresholds.keys()[0]

    if args.multiprocess:

        func = partial(structure_sizes, in_fname, args.multiprocess, var_name, ci.RESCALE, axis, ci.PERIODIC)
        params = product(t_range, var_range, thresholds.values())     
        pool = mp.Pool(processes=args.nproc)

        res = pool.map(func, params)
        pool.close()

    else:

        d = Dataset(in_fname, 'r')
        # out_fname = '%s/%s_structure_sizes_%s.pickle' % (out_dir, var_key, k)
        params = product(t_range, var_range, thresholds.values())
        res = []
        for p in params:
            res.append(structure_sizes(d, args.multiprocess, var_name, ci.RESCALE, axis, ci.PERIODIC, p))

        print ("Done.")
        d.close()

    # Reorganize data for writing
    sizes = defaultdict(dict)
    for t, v in product(t_range, var_range):
        sizes[t][v] = [x['sizes'] for x in res if x['time'] == t and x['var'] == v][0]

    out_fname = '%s/%s_structure_sizes_%s.pickle' % (out_dir, var_key, k)
    with open(out_fname, 'w') as f:
        pickle.dump(sizes, f)
