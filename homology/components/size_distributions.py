#!/usr/bin/env python

from __future__ import print_function
import numpy as np
from netCDF4 import Dataset
import argparse
import os
from homology.sampler import Sampler
# import json
import pickle
import config
from itertools import product
from functools import partial
import multiprocessing as mp

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


def structure_sizes_multiprocess(input_fname, varname, rescale, params):
    """
    Given a dataset, rescale (if necessary), split into connected components,
    and aggregate their sizes. Return as a list of integers. 
    IN
        input_fname  : path to netCDF file containing the desired data
        varname      : variable name
        params       : list of [time, height, threshold]
        **kwargs     : other parameters; rescale? If true, then it also contains 
                       offset and scale values.
    OUT
        list of integers representing the sizes of the connected components.

    """
    # Unpack parameters
    # input_fname, varname, t, z, threshold, rescale = params
    t, z, threshold = params
    # print ('Timestep %d, Height %d' % (t, z))
    # Open dataset connection
    dataset = Dataset(input_fname, 'r')

    data = dataset.variables[varname][t, z, :, :]
    if rescale:
        data = dataset.variables[varname].var_add_offset + \
               dataset.variables[varname].var_scale_factor * data

    # Close dataset connection
    dataset.close()

    if threshold == 'mean':
        sampler = Sampler(data > np.mean(data))
        sampler.connected_components()
    else:
        sampler = Sampler(data > threshold)
        sampler.connected_components()

    sizes = [len(v) for v in sampler.components.values()]
    # print sizes
    return {'time' : t,
            'height' : z,
            'threshold' : threshold,
            'sizes' : sizes}


if __name__ == '__main__':
    """
    Load file, separate components, calculate sizes.
    """
    args = process_arguments()
    print (args)

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

    # Check if filename contains variable name
    if 'VARIABLE' in ci.FILENAME:
        ci.FILENAME = ci.FILENAME.replace('VARIABLE', ci.VARIABLES[var_key])
    
    # Check if size needs to be included as well - this should only be the case for DNS dataset
    if 'SIZE' in ci.FILENAME:
        ci.FILENAME = ci.FILENAME.replace('SIZE', ci.SIZES[simulation])
    # Append simulation name as well.
    in_fname = '%s/%s/%s' % (ci.FILES, simulation, ci.FILENAME)
    
    thresholds = ci.THRESHOLDS

    print(in_fname)
    print(out_dir)

    # Open file connection to read parameters
    d = Dataset(in_fname, 'r')
    # Time range
    t_range = range(d.dimensions[ci.DIMENSIONS['t']].size)[:10]
    # print (t_range)
    # Variable range. Given by the first character in the case of velocity vector component (z_wind, etc.)
    var_direction = var_key[0]
    # print (var_direction)
    axis = d[var_name].dimensions.index(
        ci.DIMENSIONS[var_direction])
    # print (axis)
    var_range = range(d.dimensions[ci.DIMENSIONS[var_direction]].size)[:10]
    # print(var_range)
    d.close()
    
    if args.multiprocess:

        print ("Multiprocess")

        func = partial(structure_sizes_multiprocess, in_fname, var_name, ci.RESCALE)
        params = product(t_range, var_range, thresholds.values())
        
        pool = mp.Pool(processes=args.nproc)
        res = pool.map(func, params)
        pool.close()

    else:

        d = Dataset(in_fname, 'r')
        print ("inside")
        for k, threshold in thresholds.items():

            out_fname = '%s/%s_structure_sizes_%s.pickle' % (out_dir, var_key, k)
            sizes = {}

            for t in t_range:
                sizes[t] = {}
                # print (idx)
                for v in var_range:
                    # Indices for slicing
                    idx = [t, slice(None, None, None), slice(None, None, None), slice(None, None, None)]
                    idx[axis] = v
                    # separate components in two-dimensional field
                    # print (' height %d  ' % v, end='')
                    X = d.variables[var_name][idx]

                    sampler = Sampler(X > threshold)
                    sampler.connected_components(periodic=None)

                    sizes[t][v] = [x for x in sampler.uf.size if x != 0]
                    # print("t = %d, z = %d" % (t, v))
                    # print (sizes[t][v])

            print ("Done.")
        d.close()


    #########################################################################
    """
    These will depend on the dataset being used : 
    """
    #########################################################################
    # input_file = '/media/licon/1650BDF550BDDBA5/CBL3D_Data_LES/%s/wrfout_d01_2009-08-05_08-00-00' % simulation
    # input_file = '/home/licon/uni-koeln/persistence/tmp/20140717_imicro2/fielddump.w.001.nc'
    # output_dir = '/home/licon/uni-koeln/tr32/stats/size_distributions/LES/%s' % simulation
    # #########################################################################

    # dataset = Dataset(input_file, 'r')

    # varname = 'w'
    # # Time and height dimensions
    # assert(dataset.variables[varname].dimensions[0] == 'time')
    # t_range = range(dataset.variables[varname].shape[0])
    # assert(dataset.variables[varname].dimensions[1] == 'zm')
    # z_range = range(dataset.variables[varname].shape[1])


    

        # Check if output file exists
        # if os.path.isfile(out_fname):
        #     # Open as read + write
        #     with open(out_fname, 'r') as f:
        #         sizes_old = pickle.load(f)
            
        #     sizes.update(sizes_old)

        #     with open(out_fname, 'w') as f:
        #         pickle.dump(sizes, f)

        # else:
        #     # Open as write
        #     with open(out_fname, 'w') as f:
        #         pickle.dump(sizes, f)

        # dataset.close()


