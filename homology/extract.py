#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
# import scipy.io

# TODO: make this command line parameter
in_file_name = '/media/licon/1650BDF550BDDBA51/Data02/wrfout_d01_2010-04-13_12.nc'

full_dataset = Dataset(in_file_name, 'r')
for k, d in full_dataset.variables.iteritems():
    print k, d, '*' * 50
