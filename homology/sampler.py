#!/usr/bin/env python

import numpy as np
from homology.util.unionfind import UnionFind
import math

class Sampler(object):
    """
    Assume the incoming array X is two-dimensional, for now.
    """
    def __init__(self, X):
        self.uf = UnionFind()
        self.X = X
        self.components = None

    def connected_components(self, connectivity='four'):
        """
        Go through the binary array X, keeping a record of the connected components. Can use either
        four- or eight-connectivity.
        """
        if connectivity == 'four':
            it = np.nditer(self.X, flags=['multi_index'])
            for x in it:
                if x:
                    self.uf.insert_objects([it.multi_index])
                    # Check for connectivity - this is essentially as in building a merge tree, except
                    # that we don't care for actual mergers, only the end result.
                    behind = map(lambda x, y: x + y, it.multi_index, (0, -1))

                    if not any(n < 0 for n in behind) and self.X[tuple(behind)]:
                        self.uf.union(it.multi_index, tuple(behind), elder_rule=True)

                    below = map(lambda x, y: x + y, it.multi_index, (-1, 0))
                    if not any(n < 0 for n in below) and self.X[tuple(below)]:
                        self.uf.union(it.multi_index, tuple(below), elder_rule=True)

            # second pass
            it = np.nditer(self.X, flags=['multi_index'])
            for x in it:
                if x:
                    self.uf.find(it.multi_index)

            # After this we can invert the parent_pointers dict and convert to array coordinates
            self.components = {}
            for k, v in self.uf.parent_pointers.iteritems():
                keys = self.components.setdefault(self.uf.num_to_objects.get(v), [])
                keys.append(self.uf.num_to_objects.get(k))
            # Convert the lists to numpy arrays
            for k, v in self.components:
                self.components[k] = np.array(v)


        elif connectivity == 'eight':
            raise NotImplementedError
        else:
            raise ValueError('Unknown connectivity type %s.' % connectivity)

    def sample(self, points, n_samples, sample_ratio, dest=None):
        """
        points : a numpy array to sample from; rows are points in n-dimensional space.
        dest : can be either None, in which case the samples are returned, or path to an output file
        """
        n_points = points.shape[0]
        sample_size = int(math.ceil(sampling_ratio * n_points))
        total_sample = []

        if sample_size == 1:
            points_sample = points[np.random.choice(n_points)]
            total_sample.append(points_sample)
        else:
            # Perform uniform sampling without replacement
            idx = np.random.choice(points.shape[0], size=sample_size, replace=False)
            points_sample = points[idx]
            total_sample.extend(points_sample)

        total_sample = np.array(total_sample)
        



