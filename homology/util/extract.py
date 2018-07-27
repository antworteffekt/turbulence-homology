#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
from itertools import islice
from ast import literal_eval
import subprocess
import shlex
import StringIO

def array_to_cubes(X, path_to_outfile):
    """ Convert a numpy array to cube file format and write to disk.
        Convention: the rows of the array will be interpreted as points in
        Euclidean space.
    """
    with open(path_to_outfile, 'w') as out:
        out.write('; Point cloud data\n; File represents shallow cumulus clouds\ndimension 2\n')
        for row in X:
            out.write('(%d,%d)\n' % (row[0], row[1]))

def cubes_to_array(path_to_infile):
    """ Inverse of previous function: convert a cube file to numpy array.
    """
    with open(path_to_infile) as f:
        points = []
        line = f.readline()
        while line:
            # ignore comments and dimension statement
            if line[0] == ';' or line[0] == 'd':
                line = f.readline()
                continue            
            # do stuff
            points.append(literal_eval(line))
            line = f.readline()
        points = np.array(points)
        return points

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


def calculate_betti_numbers(fname, dims):
    """
    Given a file with a list of cubes, calculate its betti numbers and return them.
    """
    # Check if the file actually contains data
    lines_command_string = 'wc -l %s' % fname
    args = shlex.split(lines_command_string)
    n_lines = subprocess.check_output(args)
    n_lines = int(n_lines.split(' ')[0])
    if n_lines > 1:
        command_string = 'chomp -w %d -w %d %s' % (dims[0], dims[1], fname)
        args = shlex.split(command_string)
        betti_numbers = subprocess.check_output(args)
        return [int(x) for x in betti_numbers.split(' ')]
    else:
        return None

def persistence_intervals_to_array(filename):
    # Given a Ripser output file, create a numpy array with the corresnpoding persistence intervals per dim.
    pi_arrays = {}

    with open(filename, 'r') as f:
        pi_str = f.read()
        pi_str = pi_str.split('\n')
        # Assume the line with value range is always in position 2
        value_range = pi_str[2]
        try:
            assert value_range.startswith('value range:')
        except AssertionError as error:
            return None

        max_value = value_range.split(',')[1][:-1]
        dim_strings = [s for s in pi_str if s.startswith('persistence intervals in dim')]
        dim_idx = [pi_str.index(s) for s in dim_strings]
        # Add the index for last position in list
        dim_idx.append(-1)
        # print pi_str

        for i in range(1, len(dim_idx)):
            
            intervals = pi_str[dim_idx[i-1]:dim_idx[i]]
            identifier = intervals.pop(0).replace(':', '')

            # Infinite intervals ????
            # TODO: deal with infinite intervals without removing them, e.g. by re-adding them after 
            # conversion to array, with the maximum of the array created (not of entire range)
            intervals = [s for s in intervals if not s.endswith(', )')]

            # Remove brackets
            intervals = map(lambda x: x.replace('[', ''), intervals)
            intervals = map(lambda x: x[:-1], intervals)
            intervals = np.array([np.fromstring(x, sep=',') for x in intervals])
            
            pi_arrays[identifier] = intervals
    
    return pi_arrays


def persistence_intervals(X, ripser_bin, dim):
    
    x_str = StringIO.StringIO()
    np.savetxt(x_str, X, delimiter=',')
    args = [ripser_bin, '--format', 'point-cloud', '--dim', str(dim)]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE, 
                         stderr=subprocess.STDOUT)
    ripser_out = p.communicate(input=x_str.getvalue())[0]
    return ripser_out


def parse_output(pi_str):
    
    persistence_intervals = {}
    pi_str = pi_str.split('\n')
    value_range = pi_str[2].split(': ')[1]
    # Remove brackets
    value_range = value_range[1:-1]
    value_range = (float(value_range.split(',')[0]), float(value_range.split(',')[1]))
    persistence_intervals['value_range'] = value_range
    
    dim_strings = [s for s in pi_str if s.startswith('persistence intervals in dim')]
    dim_idx = [pi_str.index(s) for s in dim_strings]
    dim_idx.append(-1)
    for i in range(1, len(dim_idx)):

        intervals = pi_str[dim_idx[i-1]:dim_idx[i]]
        identifier = intervals.pop(0).replace(':', '')[-1]
        inf_intervals = [s for s in intervals if s.endswith(', )')]
        intervals = [s for s in intervals if not s.endswith(', )')]
        
        # Remove brackets
        intervals = map(lambda x: x.replace('[', ''), intervals)
        intervals = map(lambda x: x[:-1], intervals)
        intervals = [(float(s.split(',')[0]), float(s.split(',')[1])) for s in intervals]
        # Add infinite intervals, if present:
        if inf_intervals:
            inf_intervals = map(lambda x: x.replace('[', ''), inf_intervals)
            inf_intervals = map(lambda x: x[:-1], inf_intervals)
            inf_intervals = [(float(s.split(',')[0]), ) for s in inf_intervals]
            intervals = inf_intervals + intervals
        persistence_intervals[identifier] = intervals
    
    return persistence_intervals


