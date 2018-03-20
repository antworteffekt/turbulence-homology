#!/usr/bin/env python

import numpy as np
import time
import shlex
import subprocess
from netCDF4 import Dataset
import types


def write_cubes(X, out_fname, pos_value=1, **kwargs):
    """
    Given a two dimensional numpy array, write the corresponding list of cubes
    to disk.

    """
    # just assume that 1 is the target value for now
    X_ = np.where(X == pos_value)
    # create the string object with the appropriate data
    X_ = str(zip(X_[0], X_[1]))
    # clean up the string a bit
    X_ = X_.replace('[', '').replace(']', '')
    X_ = X_.replace('), ', ')\n')
    X_ = 'dimension 2\n' + X_
    # write to file
    outfile = open(out_fname, 'w')
    outfile.write(X_)
    outfile.close()


def write_cubes_timeblock(X, pos_value, out_fname, **kwargs):
    """
    Given a three dimensional numpy array write the corresponding list of cubes
    to disk.
    """
    # just assume that 100 is the target value for now
    X_ = np.where(X == pos_value)
    # create the string object with the appropriate data
    X_ = str(zip(X_[0], X_[1], X_[2]))
    # clean up the string a bit
    X_ = X_.replace('[', '').replace(']', '')
    X_ = X_.replace('), ', ')\n')
    X_ = 'dimension 3\n' + X_
    # write to file
    outfile = open(out_fname, 'w')
    outfile.write(X_)
    outfile.close()


def calculate_betti_numbers(fname):
    """
    Given a file with a list of cubes, calculate its betti numbers and return them.
    """
    # Check if the file actually contains data
    lines_command_string = 'wc -l %s' % fname
    args = shlex.split(lines_command_string)
    n_lines = subprocess.check_output(args)
    n_lines = int(n_lines.split(' ')[0])
    if n_lines > 1:
        command_string = 'chomp %s' % fname
        args = shlex.split(command_string)
        betti_numbers = subprocess.check_output(args)
        return [int(x) for x in betti_numbers.split(' ')]
    else:
        return None


def _binarize_array(original_data, threshold):
    """
    Given a Numpy array, either calculate a threshold value or use
    the one provided to convert it into a binary array, to be used later.
    """
    X = np.array(original_data)
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
