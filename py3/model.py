#!/usr/bin/env python

from builtins import object
import dill
import json


class HomologyModel(object):

    """
    Generate homological information from a given dataset in netCDF format.
    """

    def __init__(self, a=3, b=10):
        """
        Cosntructor.

        Parameters
        ----------

        a : int
            The spacing for the initial points of the timeblocks.
        b : int
            Length of timeblocks.
        """
        self.a = a
        self.b = b

    def save(self, file_name):
        dill.dump(self, open(file_name, 'w'))

    def load(self, file_name):
        self.__dict__.update(dill.load(open(file_name, 'r')).__dict__)

    def load_homology_data(self, file_name):
        # This considers only loading in the dictionary with the results of the
        # homological computations from a local json file.
        self.homology_data = json.load(open(file_name, 'r'))