#!/usr/bin/env python
import numpy as np
from netCDF4 import Dataset
import os
from homology.util.extract import array_to_cubes, cubes_to_array
import math
import subprocess

def sample_connected_components(cube_field, out_dir, proportional_sampling=True, sampling_proportion=0.05):
    """ 
    Given a dataset, compute its connected components and geenerate a random sample from each
    cube_field : numpy array, assumed binary
    out_dir : filesystem directory to store cube files in
    """
    cube_points = np.where(cube_field)
    cube_points = np.array((cube_points[0], cube_points[1])).T
    # If there are no cloud points, return None
    if cube_points.shape[0] == 0:
        return None
    else:
        # Check output directory
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        filelist = [f for f in os.listdir(out_dir) if f.endswith('.cub')]
        for f in filelist:
            os.remove(os.path.join(out_dir, f))
        os.chdir(out_dir)
        # Write data to disk
        out_fname = 'points.cub'
        array_to_cubes(cube_points, out_fname)
        # Call CHomP program to separate connected components
        args = ('psetconn %s component' % out_fname).split()
        # subprocess.call(args)
        subprocess.check_output(args, stderr=subprocess.STDOUT)

        # List all files produced
        filelist = [f for f in os.listdir(out_dir) if f.startswith('component')]
        total_sample = []

        # For each component, create its points array and get a sample
        for f in filelist:
            # in_fname = os.path.join(out_dir, f)
            points = cubes_to_array(f)
            n_points = points.shape[0]
            if proportional_sampling:
                # Sample a number of points proportional to component size
                sample_size = int(math.ceil(sampling_proportion * n_points))
                if sample_size == 1:
                    points_sample = points[np.random.choice(n_points)]
                    total_sample.append(points_sample)
                else:
                    # Perform uniform sampling without replacement
                    idx = np.random.choice(points.shape[0], size=sample_size, replace=False)
                    points_sample = points[idx]
                    total_sample.extend(points_sample)
            else:
                # TODO: this
                raise NotImplementedError
        # End for
        total_sample = np.array(total_sample)
        # TODO: convert points to physical coordinates? Not clear if this should be done here.
        return total_sample


def project_2d_cloud_field(dataset, timestep, variable_name):
    """
    Input: 
        dataset : netCDF4 dataset object
        timestep : integer value
    Output: 2-dimensional projection of cloud field as numpy array for given timestep, given as
    boolean array
    """
    height_dimension = 1
    # Assume the height coordinate is in position 1, otherwise projection will happen in the wrong axis
    assert dataset.variables[variable_name].dimensions[height_dimension] == 'zt'

    projected_field = np.amax(dataset.variables[variable_name][timestep, :], 0)
    projected_field = (projected_field > 0)
    return projected_field


def save_sample(dataset, sample, out_fname):
    """
    Writes both the gridpoints sampled as well as the physical coordinates.
    """
    sample_grid_fname = '%s_gridpoints.csv' % out_fname
    sample_location_fname = '%s_physicalpoints.csv' % out_fname

    # Convert grid points to physical coordinates
    yt = dataset.variables['yt'][:]
    xt = dataset.variables['xt'][:]
    # Do the map manually
    points = np.zeros(sample.shape)
    i = 0
    for row in sample:
        points[i,:] = [yt[row[0]], xt[row[1]]]
        i += 1

    with open(sample_grid_fname, 'w+') as f:
        np.savetxt(f, sample, delimiter=',')
    with open(sample_location_fname, 'w+') as f:
        np.savetxt(f, points, delimiter=',')




