#!/usr/bin/env python

import numpy as np
import gc

class HomologyModel(object):

    """Obtain, for given NetCDF data files, the corresponding homological invariants
    (Betti numbers) and represent as a time series.

    Parameters
    ----------

    a : int
        The interval size between succesive time blocks.

    b : int
        The actual time block length.




    """

    def __init__(self, a=5, b=10):
        
        self.a = a
        self.b = b


    def write_cubes(X, threshold, out_fname, **kwargs):
        """
        Given a two dimensional numpy array, write the corresponding list of cubes
        to disk.

        """
        # just assume that 100 is the target value for now
        X_ = np.where(X == 100)
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
        Given a three dimensional numpy array write the corresponding 
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