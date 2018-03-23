#!/usr/bin/env python

import numpy as np
from homology.util.unionfind import UnionFind
import math

class Sampler(object):
    """
    Assume the incoming array X is two-dimensional, for now.
    TODO: implement destructor (?)
    """
    def __init__(self, X):
        self.uf = UnionFind()
        self.X = X
        self.components = None
        self.sample_ratio = 0.05

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
            inverted_parents = {}
            for k, v in self.uf.parent_pointers.iteritems():
                keys = inverted_parents.setdefault(self.uf.num_to_objects.get(v), [])
                keys.append(self.uf.num_to_objects.get(k))
            # Convert the lists to numpy arrays
            for k, v in inverted_parents.iteritems():
                self.components[k] = np.array(v)


        elif connectivity == 'eight':
            raise NotImplementedError
        else:
            raise ValueError('Unknown connectivity type %s.' % connectivity)

    def sample(self, points, sample_ratio=0.05, dest=None):
        """
        points : a numpy array to sample from; rows are points in n-dimensional space.
        """
        n_points = points.shape[0]
        sample_size = int(math.ceil(sample_ratio * n_points))
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
        return total_sample

    def sample_components(self, n_samples=10, dest=None):
        """
        n_samples : int
        dest : if None values are returned, else it's path to output files (according to predefined format)
        """
        if self.components is None:
            self.connected_components()
        for component in self.components.values():

            for i in range(n_samples):
                sample_i = self.sample(component, self.sample_ratio)
                # print "Computed "
                # print sample_i
