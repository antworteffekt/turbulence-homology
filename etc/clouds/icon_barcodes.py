#!/usr/bin/env python
import os
import pickle
import numpy as np
import logging
from numpy.random import RandomState
from argparse import ArgumentParser
from netCDF4 import Dataset
from homology.clouds.icon import connected_components, sample_points, haversine
from scipy.spatial.distance import pdist, squareform
from homology.barcode import Barcode
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s',
                    filename='icon_barcodes.log', level=logging.INFO)


def parse_args():

    ap = ArgumentParser()
    ap.add_argument('-s', help='Simulation name.', default='controlA')
    ap.add_argument('--n-samples', help='Number of samples to compute for each\
        field.', default=10, type=int)
    ap.add_argument('-r', help='Sampling ratio.', default=0.05, type=float,
                    choices=[0.001, 0.03, 0.05])
    args = ap.parse_args()
    return args


def main(args):

    THRESHOLD = 1e-08
    prng = RandomState(17)

    if args.r == 0.001:
        ratio = '001'
    elif args.r == 0.03:
        ratio = '03'
    elif args.r == 0.05:
        ratio = '05'
    else:
        raise ValueError('Invalid sample ratio: ', args.r)

    data_dir = '/home/licon/uni-koeln/persistence/ICONdata/'
    out_root = os.path.join(data_dir, 'barcodes')
    out_dir = os.path.join(out_root, ratio, args.s)

    fname = os.path.join(data_dir, '%s_cloudwater.nc' % args.s)
    if not os.path.isfile(fname):
        raise NameError('Invalid simulation name: %s' % args.s)

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    grid_fname = os.path.join(data_dir, 'grid.nc')
    grid = Dataset(grid_fname, 'r')
    x, y = grid['clon'][:], grid['clat'][:]

    # The files have 60 timesteps.
    names = [f for f in os.listdir(out_dir) if f.endswith('pickle')]
    t_exist = [int(f.split('.')[0][1:]) for f in names]
    t_range = [t for t in range(60) if t not in t_exist]
    logging.info('Processing simulation %s, sampling ratio %s' % (args.s, ratio))
    # Save all barcodes for each timestep in a pickle file.
    for t in t_range:
        fname_out = 't%d.pickle' % t
        # Load cloud water data
        dataset = Dataset(fname, 'r')
        # Vertical maximum of cloud water
        qc = np.max(dataset['qc'][t], axis=0)
        dataset.close()
        cloud_cells = np.asarray(qc > THRESHOLD).nonzero()[0]
        if cloud_cells.size == 0:
            continue
        logging.info('%s : timestep %d' % (args.s, t))
        # Split connected components
        components = connected_components(qc > THRESHOLD, grid_fname)

        barcodes_t = []
        # Obtain n point cloud samples, and compute a barcode for each one.
        for i in range(args.n_samples):
            points = []
            for component in np.unique(components.parent):
                component_cells = cloud_cells[components.parent == component]
                sample = sample_points(prng, component_cells,
                                       sample_ratio=args.r)
                points.extend(sample)
            # X = Array with point coordinates (in radians).
            X = np.array((x[points], y[points])).T
            if X.shape[0] < 2:
                continue
            distances = squareform(pdist(X, haversine))
            barcode = Barcode()
            barcode.compute(distances, data_format='distance')
            barcodes_t.append(barcode)

        # Save barcodes for timestep t
        with open(os.path.join(out_dir, fname_out), 'wb') as f:
            pickle.dump(barcodes_t, f)

if __name__ == '__main__':
    args = parse_args()
    main(args)
