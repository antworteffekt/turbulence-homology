#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
from itertools import islice
import re

def persistence_intervals_to_array(filename):
    # Given a Ripser output file, create a numpy array with the corresnpoding persistence intervals per dim.
    pi_arrays = {}

    with open(filename, 'r') as f:
        pi_str = f.read()
        pi_str = pi_str.split('\n')
        # Assume the line with value range is always in position 2
        value_range = pi_str[2]
        assert value_range.startswith('value range:')
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
            # intervals = map(lambda x: x.replace(", )", ", %s)" % max_value), intervals)
            # TODO: deal with infinite intervals without removing them, e.g. by re-adding them after 
            # conversion to array, with the maximum of the array created (not of entire range)
            intervals = [s for s in intervals if not s.endswith(', )')]

            # Remove brackets
            intervals = map(lambda x: x.replace('[', ''), intervals)
            intervals = map(lambda x: x[:-1], intervals)
            intervals = np.array([np.fromstring(x, sep=',') for x in intervals])
            
            pi_arrays[identifier] = intervals
    
    return pi_arrays




