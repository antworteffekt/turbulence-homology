#!/usr/bin/env python
import numpy as np
from homology.util.unionfind import UnionFindOpt as UF
from math import ceil
from netCDF4 import Dataset


def connected_components(qc, grid_fname):
    """Return the union-find structure after iterating over a binary array
    and joining cells for the connected components.

    Grid file is also required."""
    grid = Dataset(grid_fname, 'r')
    n_cells = np.sum(qc)
    cell_idx = np.asarray(qc).nonzero()[0]
    uf = UF(n=n_cells)
    # maintain index to locate current cell in uf array
    i = 0
    for cell in cell_idx:
        # Shift indices to match the cell_index field
        neighbors = grid['neighbor_cell_index'][:, cell] - 1
        neighbors = neighbors[neighbors != -1]
        qc_neighbors = qc[neighbors]
        for n, qc_n in zip(neighbors, qc_neighbors):
            if qc_n:
                neighbor = np.where(cell_idx == n)[0][0]
                neighbor_parent = uf.find(neighbor)
                uf.union(i, neighbor_parent)
        i += 1
    grid.close()
    # run uf once more to compress all paths
    for i in range(n_cells):
        uf.find(i)
    return uf


def sample_points(points_array, sample_ratio=0.05):
    n_points = points_array.shape[0]
    sample_size = ceil(sample_ratio * n_points)
    sample = np.random.choice(points_array, size=sample_size, replace=False)
    return(sample)


def haversine(p1, p2):

    R = 6371.0
    lon1, lat1 = p1
    lon2, lat2 = p2

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    h = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    distance = 2 * R * np.arcsin(np.sqrt(h))
    return distance
