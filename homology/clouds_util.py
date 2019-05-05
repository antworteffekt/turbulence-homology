#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
from homology.sampler import Sampler
from homolgoy.extract import persistence_intervals, parse_output

class CloudField(object):
    """
    Cloud field object (two dimensional), with a series of samples and their corresponding VR barcodes.
    """

    def __init__(self, simulation, timestep):
        self.simulation = simulation
        self.timestep = timestep
        self.ql = None
        self.xs = None
        self.ys = None
        self.n_samples = 10
        self.samples = []
        self.barcodes = []

    def read_from_file(self, fname):
        """
        Read from projected QL field, a two-dimensional array.
        """
        dataset = Dataset(fname, 'r')
        self.ql = dataset['ql'][self.timestep,:]
        self.xs = dataset['xt'][:]
        self.ys = dataset['ys'][:]
        dataset.close()

    def sample_points(self):
        """
        Sample points from the connected components (clouds) of the ql field.
        """
        sampler = Sampler(self.ql > 0)
        sampler.connected_components(connectivity='four')
        self.samples = sampler.sample_components(n_samples=self.n_samples)

    def barcodes(self, max_dim=1):
        """
        Compute Vietoris-Rips barcodes via Ripser and store output.
        """
        # Location of Ripser binary
        ripser_bin = '/usr/local/bin/ripser/ripser'
        for s in self.samples:
            barcode = Barcode(max_dim=max_dim)
            intervals = persistence_intervals(s, ripser_bin, max_dim)
            barcode.parse_output(intervals)
            barcodes.append(barcode)

