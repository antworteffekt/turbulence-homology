#!/usr/bin/env python

import numpy as np
import time
import shlex
import subprocess
from netCDF4 import Dataset
import types


def binarize_array(X, threshold):
    """
    Given a Numpy array, either calculate a threshold value or use
    the one provided to convert it into a binary array, to be used later.
    """
    if isinstance(threshold, types.FunctionType):
        threshold_value = threshold(X)
    else:
        threshold_value = threshold
    X[X >= threshold_value] = np.max(X)
    X[X < threshold_value] = 0
    X[X == np.max(X)] = 1
    return X


def vertical_profiles(filename, variable, out_dir, z_dim_name, t_dim_name, **kwargs):
    dataset_full = Dataset(filename, 'r')
    dataset = dataset_full[variable]
    z_range = dataset.shape[dataset.dimensions.index(z_dim_name)]
    t_range = dataset.shape[dataset.dimensions.index(t_dim_name)]
    # t_range = 20

    start_time = time.time()

    profiles_0_pos = np.zeros((t_range, z_range))
    profiles_1_pos = np.zeros((t_range, z_range))
    # profiles_0_neg = np.zeros((t_range, z_range))
    # profiles_1_neg = np.zeros((t_range, z_range))

    threshold = kwargs.get('threshold', 0)
    out_fname = out_dir + 'tmp.txt'

    for t in range(t_range):
        print "Timestep %d" % t
        Xt = dataset[t, :, :, :]
        for z in range(z_range):
            X = Xt[z, :, :]
            X = _binarize_array(X, threshold=threshold)

            # "positive" domain: X == 1
            write_cubes(X, out_fname, pos_value=1)
            betti = calculate_betti_numbers(out_fname)
            if betti is not None:
                profiles_0_pos[t, z] = betti[0]
                profiles_1_pos[t, z] = betti[1]

            # "negative" domain: X == 0
            # write_cubes(X, out_fname, pos_value=0)
            # betti = calculate_betti_numbers(out_fname)
            # if betti is not None:
            #     profiles_0_neg[t, z] = betti[0]
            #     profiles_1_neg[t, z] = betti[1]

        elapsed_time = time.time() - start_time
        print "done. Time elapsed: %.3f" % elapsed_time

    # return (profiles_0_pos, profiles_1_pos, profiles_0_neg, profiles_1_neg)
    return (profiles_0_pos, profiles_1_pos)
    # return (profiles_0_neg, profiles_1_neg)


def xz_profiles(filename, variable, y_dim_name, t_dim_name, **kwargs):
    dataset_full = Dataset(filename, 'r')
    dataset = dataset_full[variable]
    y_range = dataset.shape[dataset.dimensions.index(y_dim_name)]
    t_range = dataset.shape[dataset.dimensions.index(t_dim_name)]

    profiles_0_pos = np.zeros((t_range, y_range))
    profiles_1_pos = np.zeros((t_range, y_range))
    profiles_0_neg = np.zeros((t_range, y_range))
    profiles_1_neg = np.zeros((t_range, y_range))

    threshold = kwargs.get('threshold', 0)

    for t in range(t_range):
        # for t in range(2):
        print "timestep %d" % t
        Xt = dataset[t, :, :, :]
        for y in range(y_range):
            X = Xt[:, :, y]
            X = _binarize_array(X, threshold=threshold)
            # "positive" domain: X == 1
            write_cubes(X, 'out/cubes/tmp.txt', pos_value=1)
            betti = calculate_betti_numbers('out/cubes/tmp.txt')
            if betti is not None:
                profiles_0_pos[t, y] = betti[0]
                profiles_1_pos[t, y] = betti[1]

            # "negative" domain: X == 0
            write_cubes(X, 'out/cubes/tmp.txt', pos_value=0)
            betti = calculate_betti_numbers('out/cubes/tmp.txt')
            if betti is not None:
                profiles_0_neg[t, y] = betti[0]
                profiles_1_neg[t, y] = betti[1]

    return (profiles_0_pos, profiles_1_pos, profiles_0_neg, profiles_1_neg)
