#!/usr/bin/env python
import os
import numpy as np
import time
from netCDF4 import Dataset
import types
from homology.util.extract import write_cubes, calculate_betti_numbers


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


def extract_data(model, variable, t):

    if model == 'dns':
        fname_root = '/media/licon/1650BDF550BDDBA5/CBL3D_Data_Cedrick/MEDIUM'
        if variable == 'buoy_flux':
            fname_w = 'cbl3d_it85900-91000_y_wind_subsampled_3d_x1-512_y1-512.nc'
            fname_b = 'cbl3d_it85900-91000_air_temperature_anomaly_subsampled_3d_x1-512_y1-512.nc'
            data_w = Dataset(os.path.join(fname_root, fname_w), 'r')
            data_b = Dataset(os.path.join(fname_root, fname_b), 'r')

            W = data_w['y_wind'][t] * data_w['y_wind'].var_scale_factor + data_w['y_wind'].var_add_offset
            B = data_b['air_temperature_anomaly'][t] * data_b['air_temperature_anomaly'].var_scale_factor + data_b['air_temperature_anomaly'].var_add_offset
            X = W * B
            return X

        elif variable == 'w':
            varname = 'y_wind'
            fname = 'cbl3d_it85900-91000_y_wind_subsampled_3d_x1-512_y1-512.nc'
            dataset = Dataset(os.path.join(fname_root, fname), 'r')

        elif variable == 'buoy':
            varname = 'air_temperature_anomaly'
            fname = 'cbl3d_it85900-91000_air_temperature_anomaly_subsampled_3d_x1-512_y1-512.nc'
            dataset = Dataset(os.path.join(fname_root, fname), 'r')

        X = dataset[varname][t] * dataset[varname].var_scale_factor + dataset[varname].var_add_offset
        dataset.close()
        return X



def vertical_profiles(filename, variable, out_fname, z_dim_name, t_value,
                      domain='pos', rescale=False, center=False, threshold=0,
                      periodic=False):
    dataset = Dataset(filename, 'r')
    # dataset = dataset_full[variable]

    z_range = dataset[variable].shape[dataset[variable].dimensions.index(z_dim_name)]

    start_time = time.time()

    profiles_0 = np.zeros((1, z_range))
    profiles_1 = np.zeros((1, z_range))

    t = t_value
    # for t in t_values:
    print("Timestep %d" % t)
    Xt = dataset[variable][t, :, :, :]
    if rescale:
        Xt = Xt * dataset[variable].var_scale_factor + dataset[variable].var_add_offset
    for z in range(z_range):
        X = np.array(Xt[z, :, :])
        if center:
            X = X - np.mean(X)
        if domain == 'pos':
            X = X > threshold
        elif domain == 'neg':
            X = X < threshold
        else:
            raise ValueError('Unknown domain type %s' % domain)
        dims = X.shape
        write_cubes(X, out_fname, pos_value=1)
        betti = calculate_betti_numbers(out_fname, periodic, dims)
        if betti is not None:
            profiles_0[0, z] = betti[0]
            profiles_1[0, z] = betti[1]

    elapsed_time = time.time() - start_time
    print("done. Time elapsed: %.3f" % elapsed_time)
    dataset.close()
    return (profiles_0, profiles_1)


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
        print( "timestep %d" % t)
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
