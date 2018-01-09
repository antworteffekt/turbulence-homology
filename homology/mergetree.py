#!/usr/bin/env python

import numpy as np
from homology.util.unionfind import UnionFind
import matplotlib.pyplot as plt
from collections import defaultdict

import pprint


class MergeTree(object):

    def __init__(self):
        self.uf = UnionFind()
        self.levels_components = {}
        self.components_level = {}
        self.mergers = {}
        self.dimension = 0

    def build_merge_tree(self, X, iter_dim=0, flip=False, **kwargs):
        """
        Given a two-dimensional array with binarized data, builds a union-find
        data structure incrementally by iterating over the specified dimension
        (physical "height"). For each array level, this keeps track of the number
        of connected components of TRUE values.
        """

        # Assume X is not binarized yet, do it in here:
        tol = kwargs.get('tol', 0.01)
        X_bin = X > tol
        # X_bin = X
        # print "X binarized"

        it = np.nditer(X_bin, flags=['multi_index'])
        self.dimension = 2
        below_joins = {}
        while not it.finished:
            current_level = it.multi_index[0]
            # print "Current level: %d" % current_level
            if it[0]:
                # print it.multi_index
                # print X[it.multi_index]
                # self.uf.insert_objects([it.multi_index])

                # Special case: level 0. Here we must assign roots manually, to
                # accurately determine the location of the maximum values.
                # if current_level == 0:
                #     print self.uf.parent_pointers
                #     print "Level 0"
                #     cells = defaultdict(list)
                #     for k, v in self.uf.parent_pointers.iteritems():
                #         cells[v].append(k)
                #     # tmp_cmp = {}
                #     print cells
                        # pass
                
                # Pass both the array coordinates and the array value to find
                self.uf.insert_objects([it.multi_index], value=X[it.multi_index])

                # Check for connectivity behind and below only - this works because
                # we assume a two-dimensional array structure
                behind = map(lambda x, y: x + y, it.multi_index, (0, -1))
                if not any(n < 0 for n in behind) and X_bin[tuple(behind)]:
                    self.uf.union(it.multi_index, tuple(behind), elder_rule=False)

                below = map(lambda x, y: x + y, it.multi_index, (-1, 0))
                if not any(n < 0 for n in below) and X_bin[tuple(below)]:
                    # Don't do the union here, just keep track of the involved
                    # elements instead.
                    # self.uf.union(it.multi_index, tuple(below), elder_rule=False)
                    # Get cell's parent in this level
                    cell_parent = self.uf.find(it.multi_index)
                    # Get root of cell below
                    below_parent = self.uf.find(tuple(below))
                    # print "*** Below tuple: %s" % str(tuple(below))
                    # print "*** Below parent: %s" % str(below_parent)
                    if cell_parent not in below_joins:
                        below_joins[cell_parent] = [tuple(below_parent)]
                    else:
                        # print "appending %s" % str(below_parent)
                        below_joins[cell_parent].append(below_parent)

            it.iternext()
            if it.finished or it.multi_index[0] != current_level:
                # print "Level %d finished, processing with weights:" % current_level
                # print self.uf.num_weights
                # print "%d components" % len(self.uf.num_weights)
                # print self.uf.parent_pointers
                # print "Level %d" % current_level
                # print "Joins to be made" 
                # pprint.pprint(below_joins)
                for level_parent, below_parent in below_joins.iteritems():
                    # print k, d, len(set(d))
                    # print level_parent
                    # print below_parent
                    if len(set(below_parent)) > 1:
                        # Get the roots for the objects below
                        # print "*** Joining group rooted at %s ***" % str(level_parent)
                        roots_below = list(set(below_parent))
                        # Use the find method on the objects first, in case one of the 
                        # roots has already been merged before.
                        id_roots_below = [self.uf.find(x) for x in roots_below]
                        id_roots_below = [self.uf.objects_to_num[x] for x in id_roots_below]
                        # print "Roots below: %s" % str(list(set(below_parent)))
                        # print "Ids : %s" % str(id_roots_below)
                        # 0th element is at once previous and maximum
                        # print "Weights before union:"
                        # print self.uf.num_weights
                        max_root_id = id_roots_below[0]
                        # print "Max root parent: %s" % str(self.uf.find(roots_below[0]))
                        # print "Value at position 0, current maximum: %f" % self.uf.num_weights[roots_below[0]]
                        iter_roots = iter(id_roots_below)
                        next(iter_roots)
                        for id_root in iter_roots:
                            # print "Value at next position: %f" % self.uf.num_weights[root]
                            # print "Prev: %s" % str(prev)
                            # print "Curr: %s" % str(root)
                            # print "Current root: %s" % str(self.uf.num_to_objects[id_root])
                            # print "Current root parent: %s" % str(self.uf.find(self.uf.num_to_objects[id_root]))
                            # print "max root id: %d" % max_root_id
                            # print "value: %f" % self.uf.num_weights[max_root_id]
                            # print "id root --- %d" % id_root
                            # print "value: %f" % self.uf.num_weights[id_root]
                            if self.uf.num_weights[max_root_id] < self.uf.num_weights[id_root]:
                                max_root_id = id_root

                        # print "*** Several roots below, maximum was selected"
                        # print "Level parent before union: %s" % str(level_parent)
                        # print "Parent below: %s" % str(self.uf.num_to_objects[max_root])
                        # Join current level group with the maximally valued root
                        max_root = self.uf.num_to_objects[max_root_id]
                        # print "Root with maximum value found --- %d : %s" % (max_root_id, str(max_root))

                        # print "Maximum value: %f" % self.uf.num_weights[max_root_id]
                        self.uf.union(level_parent, max_root, elder_rule=True)
                        # print "Weights after union:"
                        # print self.uf.num_weights

                        # Join the other non-maximally valued roots:
                        # print "Below parents: %s" % str(roots_below)
                        missing_roots = [root for root in roots_below if root != max_root]
                        # print "Missing roots: %s" % str(missing_roots)
                        for root in missing_roots:
                            # print "Joining root %s : %d with maximum %s : %d" % (str(root), self.uf.objects_to_num[root], str(max_root), max_root_id)
                            self.uf.union(max_root, root, elder_rule=True)
                        # print "Weights after union of remaining roots:"
                        # print self.uf.num_weights
                        # print "\n"
                        # print "After union ... %s" % str(self.uf.find(level_parent))
                    else:  # there is just one root below
                        # print "*** There is just one root below"
                        # print "Level parent before union: %s" % str(level_parent)
                        # print "Parent below: %s" % str(below_parent[0])
                        # print self.uf
                        self.uf.union(level_parent, below_parent[0], elder_rule=True)
                        # print "After union ... %s" % str(self.uf.find(level_parent))
                        # pass

                    # print "Weights after merging, etc. (%d components):" % len(self.uf.num_weights) 
                    # print self.uf.num_weights

                level_coords = [x for x in self.uf.objects_to_num.keys() if x[0] == current_level]
                self.levels_components[current_level] = set([self.uf.find(x) for x in level_coords])
                for key in self.levels_components[current_level]:
                    if key in self.components_level:
                        self.components_level[key].append(current_level)
                    else:
                        self.components_level[key] = [current_level]

                # Check if any of the components from previous level "disappeared"
                if flip:
                    prev_level = min(X_bin.shape[iter_dim] - 1, current_level + 1)
                else:
                    prev_level = max(0, current_level - 1)
                # print "Previous level: %d\n" % prev_level
                level_difference = self.levels_components[prev_level] - self.levels_components[current_level]
                if len(level_difference) > 0:
                    for cell in level_difference:
                        new_parent = self.uf.find(cell)
                        if new_parent != cell:
                            self.mergers[cell] = (prev_level, new_parent)
                # Reset the below_joins dictionary
                below_joins.clear()

    def build_merge_tree_3d(self, X, iter_dim=0, flip=False, **kwargs):
        """
        This is ridiculous - find a way to join this with the above function
        """
        # Assume X is not binarized yet, do it in here:
        tol = kwargs.get('tol', 0.01)
        X_bin = X > tol
        # X_bin = X
        # print "X binarized"

        it = np.nditer(X_bin, flags=['multi_index'])
        self.dimension = 3
        below_joins = {}
        while not it.finished:
            current_level = it.multi_index[0]
            # print "Current level: %d" % current_level
            if it[0]:
                # print it.multi_index
                # print X[it.multi_index]
                # self.uf.insert_objects([it.multi_index])

                # Special case: level 0. Here we must assign roots manually, to
                # accurately determine the location of the maximum values.
                # if current_level == 0:
                #     print self.uf.parent_pointers
                #     print "Level 0"
                #     cells = defaultdict(list)
                #     for k, v in self.uf.parent_pointers.iteritems():
                #         cells[v].append(k)
                #     # tmp_cmp = {}
                #     print cells
                        # pass
                
                # Pass both the array coordinates and the array value to find
                self.uf.insert_objects([it.multi_index], value=X[it.multi_index])

                # Check for connectivity behind and below only - this works because
                # we assume a two-dimensional array structure
                behind = map(lambda x, y: x + y, it.multi_index, (0, 0, -1))
                if not any(n < 0 for n in behind) and X_bin[tuple(behind)]:
                    self.uf.union(it.multi_index, tuple(behind), elder_rule=False)

                side = map(lambda a, b: a + b, it.multi_index, (0, -1, 0))
                if not any(n < 0 for n in side) and X[tuple(side)]:
                    self.uf.union(it.multi_index, tuple(side), elder_rule=False)

                below = map(lambda x, y: x + y, it.multi_index, (-1, 0, 0))
                if not any(n < 0 for n in below) and X_bin[tuple(below)]:
                    # Don't do the union here, just keep track of the involved
                    # elements instead.
                    # self.uf.union(it.multi_index, tuple(below), elder_rule=False)
                    # Get cell's parent in this level
                    cell_parent = self.uf.find(it.multi_index)
                    # Get root of cell below
                    below_parent = self.uf.find(tuple(below))
                    # print "*** Below tuple: %s" % str(tuple(below))
                    # print "*** Below parent: %s" % str(below_parent)
                    if cell_parent not in below_joins:
                        below_joins[cell_parent] = [tuple(below_parent)]
                    else:
                        # print "appending %s" % str(below_parent)
                        below_joins[cell_parent].append(below_parent)

            it.iternext()
            if it.finished or it.multi_index[0] != current_level:
                # print "Level %d finished, processing with weights:" % current_level
                # print self.uf.num_weights
                # print "%d components" % len(self.uf.num_weights)
                # print self.uf.parent_pointers
                # print "Level %d" % current_level
                # print "Joins to be made" 
                # pprint.pprint(below_joins)
                for level_parent, below_parent in below_joins.iteritems():
                    # print k, d, len(set(d))
                    # print level_parent
                    # print below_parent
                    if len(set(below_parent)) > 1:
                        # Get the roots for the objects below
                        # print "*** Joining group rooted at %s ***" % str(level_parent)
                        roots_below = list(set(below_parent))
                        # Use the find method on the objects first, in case one of the 
                        # roots has already been merged before.
                        id_roots_below = [self.uf.find(x) for x in roots_below]
                        id_roots_below = [self.uf.objects_to_num[x] for x in id_roots_below]
                        # print "Roots below: %s" % str(list(set(below_parent)))
                        # print "Ids : %s" % str(id_roots_below)
                        # 0th element is at once previous and maximum
                        # print "Weights before union:"
                        # print self.uf.num_weights
                        max_root_id = id_roots_below[0]
                        # print "Max root parent: %s" % str(self.uf.find(roots_below[0]))
                        # print "Value at position 0, current maximum: %f" % self.uf.num_weights[roots_below[0]]
                        iter_roots = iter(id_roots_below)
                        next(iter_roots)
                        for id_root in iter_roots:
                            # print "Value at next position: %f" % self.uf.num_weights[root]
                            # print "Prev: %s" % str(prev)
                            # print "Curr: %s" % str(root)
                            # print "Current root: %s" % str(self.uf.num_to_objects[id_root])
                            # print "Current root parent: %s" % str(self.uf.find(self.uf.num_to_objects[id_root]))
                            # print "max root id: %d" % max_root_id
                            # print "value: %f" % self.uf.num_weights[max_root_id]
                            # print "id root --- %d" % id_root
                            # print "value: %f" % self.uf.num_weights[id_root]
                            if self.uf.num_weights[max_root_id] < self.uf.num_weights[id_root]:
                                max_root_id = id_root

                        # print "*** Several roots below, maximum was selected"
                        # print "Level parent before union: %s" % str(level_parent)
                        # print "Parent below: %s" % str(self.uf.num_to_objects[max_root])
                        # Join current level group with the maximally valued root
                        max_root = self.uf.num_to_objects[max_root_id]
                        # print "Root with maximum value found --- %d : %s" % (max_root_id, str(max_root))

                        # print "Maximum value: %f" % self.uf.num_weights[max_root_id]
                        self.uf.union(level_parent, max_root, elder_rule=True)
                        # print "Weights after union:"
                        # print self.uf.num_weights

                        # Join the other non-maximally valued roots:
                        # print "Below parents: %s" % str(roots_below)
                        missing_roots = [root for root in roots_below if root != max_root]
                        # print "Missing roots: %s" % str(missing_roots)
                        for root in missing_roots:
                            # print "Joining root %s : %d with maximum %s : %d" % (str(root), self.uf.objects_to_num[root], str(max_root), max_root_id)
                            self.uf.union(max_root, root, elder_rule=True)
                        # print "Weights after union of remaining roots:"
                        # print self.uf.num_weights
                        # print "\n"
                        # print "After union ... %s" % str(self.uf.find(level_parent))
                    else:  # there is just one root below
                        # print "*** There is just one root below"
                        # print "Level parent before union: %s" % str(level_parent)
                        # print "Parent below: %s" % str(below_parent[0])
                        # print self.uf
                        self.uf.union(level_parent, below_parent[0], elder_rule=True)
                        # print "After union ... %s" % str(self.uf.find(level_parent))
                        # pass

                    # print "Weights after merging, etc. (%d components):" % len(self.uf.num_weights) 
                    # print self.uf.num_weights

                level_coords = [x for x in self.uf.objects_to_num.keys() if x[0] == current_level]
                self.levels_components[current_level] = set([self.uf.find(x) for x in level_coords])
                for key in self.levels_components[current_level]:
                    if key in self.components_level:
                        self.components_level[key].append(current_level)
                    else:
                        self.components_level[key] = [current_level]

                # Check if any of the components from previous level "disappeared"
                if flip:
                    prev_level = min(X_bin.shape[iter_dim] - 1, current_level + 1)
                else:
                    prev_level = max(0, current_level - 1)
                # print "Previous level: %d\n" % prev_level
                level_difference = self.levels_components[prev_level] - self.levels_components[current_level]
                if len(level_difference) > 0:
                    for cell in level_difference:
                        new_parent = self.uf.find(cell)
                        if new_parent != cell:
                            self.mergers[cell] = (prev_level, new_parent)
                # Reset the below_joins dictionary
                below_joins.clear()

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
        # ax.invert_xaxis()
        fig.show()

    # def _binarize_array(self, X, tol=0.01):
