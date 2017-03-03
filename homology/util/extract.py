#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
# import scipy.io

# TODO: make this command line parameter
in_file_name = '/media/licon/1650BDF550BDDBA52/Data02/wrfout_d01_2010-04-13_12.nc'
OUT_DIR = 'data/data02/'

full_dataset = Dataset(in_file_name, 'r')
# for k, d in full_dataset.variables.iteritems():
#     print k, d, '*' * 50 
variables = {
    'T': 'temperature',
    'QVAPOR': 'qvapor',
    'Q2': 'q2',
    'U': 'u',
    'V': 'v',
    'W': 'w'
}

# choose variables of interest - here I assume they contain the full 4 dimensions
temperature = atmosphere.variables['T'][:,:,:,:]
qvapor = atmosphere.variables['QVAPOR'][:,:,:,:]
q2 = atmosphere.variables['Q2'][:,:,:] # humidity at height 2 exactly
u = atmosphere.variables['U'][:,:,:,:] # x-component of wind
v = atmosphere.variables['V'][:,:,:,:] # y-component of wind
w = atmosphere.variables['W'][:,:,:,:] # z-component of wind

# choose variables of interest - here I assume they contain the full 4 dimensions
temperature = atmosphere.variables['T'][:,:,:,:]
qvapor = atmosphere.variables['QVAPOR'][:,:,:,:]
q2 = atmosphere.variables['Q2'][:,:,:] # humidity at height 2 exactly
u = atmosphere.variables['U'][:,:,:,:] # x-component of wind
v = atmosphere.variables['V'][:,:,:,:] # y-component of wind
w = atmosphere.variables['W'][:,:,:,:] # z-component of wind