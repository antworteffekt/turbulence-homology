#!/usr/bin/env python
import os
import numpy as np
from argparse import ArgumentParser



def parse_args():

    ap = ArgumentParser()
    ap.add_argument('-s', help='Simulation name.', default='controlA')
    ap.add_argument('--n-samples', help='Number of samples to compute for each\
        field.', default=10, type=int)
    ap.add_argument('-r', help='Sampling ratio.', default=0.05, type=float)
    args = ap.parse_args()
    return args


def main():

    args = parse_args()

    print(args)

if __name__ == '__main__':
    main()
