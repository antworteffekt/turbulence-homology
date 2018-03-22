#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import argparse

description = """
Given a dataset, generate point samples which will be used to build Vietoris-Rips filtrations.
"""

CLI = argparse.ArgumentParser(
    description='Generate point samples to be used as point cloud input for persistent homology')
CLI.add_argument(
    'data',
    help='Dataset to use (NetCDF file)')
CLI.add_argument(
    'varname',
    help='Variable name to use for sampling')

args = CLI.parse_args()

def main():
    pass

main()