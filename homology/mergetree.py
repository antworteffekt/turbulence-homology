#!/usr/bin/env python

import numpy as np
from homology.util.unionfind import UnionFind
import matplotlib.pyplot as plt


class MergeTree(object):

    def __init__(self):
        self.uf = UnionFind()
        self.levels_components = {}
        self.components_level = {}
        self.mergers = {}
        self.dimension = 0

    def build_merge_tree(self, X, iter_dim=0, flip=False):
        """
        Given a two-dimensional array with binarized data, builds a union-find
        data structure incrementally by iterating over the specified dimension
        (physical "height"). For each array level, this keeps track of the number
        of connected components of TRUE values.
        """
        it = np.nditer(X, flags=['multi_index'])
        self.dimension = 2
        while not it.finished:
            current_level = it.multi_index[0]

            if it[0]:
                # print it.multi_index
                self.uf.insert_objects([it.multi_index])
                # Check for connectivity behind and below only - this works because
                # we assume a two-dimensional array structure
                behind = map(lambda x, y: x + y, it.multi_index, (0, -1))
                if not any(n < 0 for n in behind) and X[tuple(behind)]:
                    self.uf.union(it.multi_index, tuple(behind), elder_rule=False)

                below = map(lambda x, y: x + y, it.multi_index, (-1, 0))
                if not any(n < 0 for n in below) and X[tuple(below)]:
                    self.uf.union(it.multi_index, tuple(below), elder_rule=False)

            it.iternext()
            if it.finished or it.multi_index[0] != current_level:
                # print "Current level: %d" % current_level
                level_coords = [x for x in self.uf.objects_to_num.keys() if x[0] == current_level]
                self.levels_components[current_level] = set([self.uf.find(x) for x in level_coords])
                for key in self.levels_components[current_level]:
                    if key in self.components_level:
                        self.components_level[key].append(current_level)
                    else:
                        self.components_level[key] = [current_level]

                # Check if any of the components from previous level "disappeared"
                if flip:
                    prev_level = min(X.shape[iter_dim] - 1, current_level + 1)
                else:
                    prev_level = max(0, current_level - 1)
                # print "Previous level: %d\n" % prev_level
                level_difference = self.levels_components[prev_level] - self.levels_components[current_level]
                if len(level_difference) > 0:
                    for cell in level_difference:
                        new_parent = self.uf.find(cell)
                        if new_parent != cell:
                            self.mergers[cell] = (prev_level, new_parent)

    def build_merge_tree_3d(self, X, iter_dim=0):
        """
        This is ridiculous - find a way to join this with the above function
        """
        it = np.nditer(X, flags=['multi_index'])
        self.dimension = 3
        while not it.finished:
            current_level = it.multi_index[0]
            if it[0]:
                self.uf.insert_objects([it.multi_index])
                behind = map(lambda x, y: x + y, it.multi_index, (0, 0, -1))
                if not any(n < 0 for n in behind) and X[tuple(behind)]:
                    self.uf.union(it.multi_index, tuple(behind), elder_rule=False)

                side = map(lambda a, b: a + b, it.multi_index, (0, -1, 0))
                if not any(n < 0 for n in side) and X[tuple(side)]:
                    self.uf.union(it.multi_index, tuple(side), elder_rule=False)

                below = map(lambda a, b: a + b, it.multi_index, (-1, 0, 0))
                if not any(n < 0 for n in below) and X[tuple(below)]:
                    self.uf.union(it.multi_index, tuple(below), elder_rule=False)

            it.iternext()
            if it.finished or it.multi_index[0] != current_level:
                # perform level-wise checking of components, etc.
                print "checking components for level %d" % current_level
                level_coords = [x for x in self.uf.objects_to_num.keys() if x[0] == current_level]
                self.levels_components[current_level] = set([self.uf.find(x) for x in level_coords])
                for key in self.levels_components[current_level]:
                    if key in self.components_level:
                        self.components_level[key].append(current_level)
                    else:
                        self.components_level[key] = [current_level]

                # check if there are any components that "disappeared"
                prev_level = max(0, current_level - 1)
                level_difference = self.levels_components[prev_level] - self.levels_components[current_level]
                if len(level_difference) > 0:
                    for cell in level_difference:
                        new_parent = self.uf.find(cell)
                        if new_parent != cell:
                            self.mergers[cell] = (prev_level, new_parent)

    def plot(self):
        fig = plt.figure()
        ax = plt.axes()
        for k, d in self.components_level.iteritems():
            ax.vlines(k[1], d[0], d[-1])
        for k, d in self.mergers.iteritems():
            xmin = min(k[1], d[1][1])
            xmax = max(k[1], d[1][1])
            ax.hlines(d[0] + 1, xmin, xmax)
            ax.vlines(k[1], d[0], d[0] + 1)
            ax.plot(d[1][1], d[0] + 1, 'r.')
        fig.show()

    # def _binarize_array(self, X, tol=0.01):
