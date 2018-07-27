#!/usr/bin/env python

import numpy as np
from homology.util.unionfind import UnionFind
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import networkx as nx
import pprint


class MergeTree(object):

    def __init__(self, save_values=False):
        self.uf = UnionFind()
        self.levels_components = {}
        self.components_level = {}
        self.mergers = {}
        self.dimension = 0
        if save_values:
            self.component_sum = {}
        else:
            self.component_sum = None

    def build_merge_tree(self, X, periodic=False, return_graph=False, 
                         iter_dim=0, flip=False, **kwargs):
        """
        Given a two-dimensional array with binarized data, builds a union-find
        data structure incrementally by iterating over the specified dimension
        (physical "height"). For each array level, this keeps track of the number
        of connected components of TRUE values.

        This uses 4-connectivity.
        TODO: implement the method for 8-connectivity.
        """

        # Assume X is not binarized yet, do it in here:

        it = np.nditer(X, flags=['multi_index'])
        self.dimension = 2
        below_joins = {}
        if return_graph:
            G = nx.DiGraph()

        while not it.finished:
            current_level = it.multi_index[0]
            if it[0]:
                # Pass both the array coordinates and the array value to find
                cell_parent = self.uf.insert_objects([it.multi_index], value=X[it.multi_index])
                print type(cell_parent)
                print 'Cell: %s  --- parent: %s' % (it.multi_index, cell_parent)
                # Check for connectivity behind and below only - this works because
                # we assume a two-dimensional array structure
                # Periodicity should affect only the lookup of the cell "behind".
                if periodic:
                    behind = map(lambda x, y: x + y, it.multi_index, (0, -1))
                    # Map the column index to its residue mod n_columns.
                    behind = (behind[0], behind[1] % X.shape[1])
                    if X[tuple(behind)]:
                        self.uf.union(it.multi_index, tuple(behind), elder_rule=False)
                else:
                    behind = map(lambda x, y: x + y, it.multi_index, (0, -1))
                    if not any(n < 0 for n in behind) and X[tuple(behind)]:
                        self.uf.union(it.multi_index, tuple(behind), elder_rule=False)

                below = map(lambda x, y: x + y, it.multi_index, (-1, 0))
                if not any(n < 0 for n in below) and X[tuple(below)]:
                    # Don't do the union here, just keep track of the involved
                    # elements instead.
                    # Get cell's parent in this level
                    cell_parent = self.uf.find(it.multi_index)
                    # Get root of cell below
                    below_parent = self.uf.find(tuple(below))
                    if cell_parent not in below_joins:
                        below_joins[cell_parent] = [tuple(below_parent)]
                    else:
                        below_joins[cell_parent].append(below_parent)

            it.iternext()
            if it.finished or it.multi_index[0] != current_level:
                for level_parent, below_parent in below_joins.iteritems():
                    if len(set(below_parent)) > 1:
                        # Get the roots for the objects below
                        roots_below = list(set(below_parent))
                        # Use the find method on the objects first, in case one of the 
                        # roots has already been merged before.
                        id_roots_below = [self.uf.find(x) for x in roots_below]
                        id_roots_below = [self.uf.objects_to_num[x] for x in id_roots_below]
                        # 0th element is at once previous and maximum
                        max_root_id = id_roots_below[0]
                        iter_roots = iter(id_roots_below)
                        next(iter_roots)
                        for id_root in iter_roots:
                            if self.uf.num_weights[max_root_id] < self.uf.num_weights[id_root]:
                                max_root_id = id_root
                        # Join current level group with the maximally valued root
                        max_root = self.uf.num_to_objects[max_root_id]
                        self.uf.union(level_parent, max_root, elder_rule=True)
                        # Join the other non-maximally valued roots:
                        missing_roots = [root for root in roots_below if root != max_root]
                        for root in missing_roots:
                            self.uf.union(max_root, root, elder_rule=True)
                    else:  # there is just one root below
                        self.uf.union(level_parent, below_parent[0], elder_rule=True)

                level_coords = [x for x in self.uf.objects_to_num.keys() if x[0] == current_level]
                level_parents = [self.uf.find(x) for x in level_coords]

                parents_points = {x : [] for x in set(level_parents)}
                for p in level_coords:
                    parents_points[self.uf.find(p)].append(p)

                self.levels_components[current_level] = set(parents_points.keys())

                for key in self.levels_components[current_level]:
                    if key in self.components_level:
                        self.components_level[key].append(current_level)
                    else:
                        # Component has not been seen before.
                        self.components_level[key] = [current_level]
                        # Add to the graph as root node if required
                        if return_graph:
                            G.add_node(key, type='root', level=current_level, points=parents_points[key])
                
                # Check if any of the components from previous level "disappeared", this signals a merger.
                if flip:
                    prev_level = min(X.shape[iter_dim] - 1, current_level + 1)
                else:
                    prev_level = max(0, current_level - 1)

                level_difference = self.levels_components[prev_level] - self.levels_components[current_level]
                if len(level_difference) > 0:
                    for cell in level_difference:
                        new_parent = self.uf.find(cell)
                        if new_parent != cell:
                            # Merger event
                            self.mergers[cell] = (prev_level, new_parent)
                            if return_graph:
                                merge_node = (prev_level + 1, new_parent[1])
                                G.add_node(merge_node, type='merge', level=current_level)
                                G.add_edge(cell, merge_node)
                                G.add_edge(new_parent, merge_node)
                        else:
                            # No merger, but add terminal node
                            end_node = (current_level, cell[1])
                            if return_graph:
                                G.add_node(end_node, type='leaf', level=current_level)
                                G.add_edge(cell, end_node)
                # Reset the below_joins dictionary
                below_joins.clear()
                # print "Level %d done" % current_level
                # print "Components present: ", self.levels_components[current_level]
            # End if
        # End while (going through horizontal levels)

        if return_graph:
            return G

    def build_merge_tree_3d(self, X, periodic=False, return_graph=False, iter_dim=0, flip=False, **kwargs):
        """
        This is ridiculous - find a way to join this with the above function
        """
        # Assume X is not binarized yet, do it in here:
        # tol = kwargs.get('tol', 0.01)
        # X_bin = X > tol

        it = np.nditer(X, flags=['multi_index'])
        self.dimension = 3
        below_joins = {}
        if return_graph:
            G = nx.DiGraph()

        while not it.finished:
            current_level = it.multi_index[0]
            if it[0]:
                # Pass both the array coordinates and the array value to find
                self.uf.insert_objects([it.multi_index], value=X[it.multi_index])

                if periodic:
                    behind = map(lambda x, y: x + y, it.multi_index, (0, 0, -1))
                    behind = (behind[0], behind[1], behind[2] % X.shape[2])
                    if X[tuple(behind)]:
                        self.uf.union(it.multi_index, tuple(behind), elder_rule=False)

                    side = map(lambda a, b: a + b, it.multi_index, (0, -1, 0))
                    side = (side[0], side[1] % X.shape[1], side[2])
                    if X[tuple(side)]:
                        self.uf.union(it.multi_index, tuple(side), elder_rule=False)
                else:
                    # Check for connectivity behind, below, and to one side - this assumes 
                    # six connectivity
                    behind = map(lambda x, y: x + y, it.multi_index, (0, 0, -1))
                    if not any(n < 0 for n in behind) and X[tuple(behind)]:
                        self.uf.union(it.multi_index, tuple(behind), elder_rule=False)

                    side = map(lambda a, b: a + b, it.multi_index, (0, -1, 0))
                    if not any(n < 0 for n in side) and X[tuple(side)]:
                        self.uf.union(it.multi_index, tuple(side), elder_rule=False)

                below = map(lambda x, y: x + y, it.multi_index, (-1, 0, 0))
                if not any(n < 0 for n in below) and X[tuple(below)]:
                    # Don't do the union here, just keep track of the involved elements instead.
                    # Get cell's parent in this level
                    cell_parent = self.uf.find(it.multi_index)
                    # Get root of cell below
                    below_parent = self.uf.find(tuple(below))
                    if cell_parent not in below_joins:
                        below_joins[cell_parent] = [tuple(below_parent)]
                    else:
                        below_joins[cell_parent].append(below_parent)

            it.iternext()
            if it.finished or it.multi_index[0] != current_level:
                # Level empty, store merger information
                for level_parent, below_parent in below_joins.iteritems():
                    if len(set(below_parent)) > 1:
                        # Get the roots for the objects below
                        roots_below = list(set(below_parent))
                        # Use the find method on the objects first, in case one of the
                        # roots has already been merged before.
                        id_roots_below = [self.uf.find(x) for x in roots_below]
                        id_roots_below = [self.uf.objects_to_num[x] for x in id_roots_below]
                        # 0th element is at once previous and maximum
                        max_root_id = id_roots_below[0]
                        iter_roots = iter(id_roots_below)
                        next(iter_roots)
                        for id_root in iter_roots:
                            if self.uf.num_weights[max_root_id] < self.uf.num_weights[id_root]:
                                max_root_id = id_root
                        # Join current level group with the maximally valued root
                        max_root = self.uf.num_to_objects[max_root_id]
                        self.uf.union(level_parent, max_root, elder_rule=True)
                        # Join the other non-maximally valued roots:
                        missing_roots = [root for root in roots_below if root != max_root]
                        for root in missing_roots:
                            self.uf.union(max_root, root, elder_rule=True)
                    else:  # there is just one root below
                        self.uf.union(level_parent, below_parent[0], elder_rule=True)

                level_coords = [x for x in self.uf.objects_to_num.keys() if x[0] == current_level]
                level_parents = [self.uf.find(x) for x in level_coords]

                parents_points = {x : [] for x in set(level_parents)}
                for p in level_coords:
                    parents_points[self.uf.find(p)].append(p)

                self.levels_components[current_level] = set(parents_points.keys())

                for key in self.levels_components[current_level]:
                    if key in self.components_level:
                        self.components_level[key].append(current_level)
                    else:
                    # Component has not been seen before: birth event => add node to graph
                        self.components_level[key] = [current_level]
                        if return_graph:
                            G.add_node(key, type='root', level=current_level, points=parents_points[key])

                # Check if any of the components from previous level "disappeared"
                if flip:
                    prev_level = min(X.shape[iter_dim] - 1, current_level + 1)
                else:
                    prev_level = max(0, current_level - 1)
                level_difference = self.levels_components[prev_level] - self.levels_components[current_level]
                if len(level_difference) > 0:
                    for cell in level_difference:
                        new_parent = self.uf.find(cell)
                        if new_parent != cell:
                            # Merger event, add node
                            self.mergers[cell] = (prev_level, new_parent)
                            if return_graph:
                                merge_node = (current_level, new_parent[1], new_parent[2])
                                G.add_node(merge_node, type='merge', level=current_level)
                                G.add_edge(cell, merge_node)
                                G.add_edge(new_parent, merge_node)
                        else:
                            # No merger, but add terminal node
                            end_node = (current_level, cell[1], cell[2])
                            if return_graph:
                                G.add_node(end_node, type='leaf', level=current_level)
                                G.add_edge(cell, end_node)
                # Reset the below_joins dictionary
                below_joins.clear()

        if return_graph:
            return G

    def structure_size(self):
        """
        Count the number of cells in each connected structure.
        Returns {parent : #(cells)}
        """
        # First invert the dict parent_pointers
        child_nodes = {}
        for k, v in self.uf.parent_pointers.iteritems():
            keys = child_nodes.setdefault(v, [])
            keys.append(k)
        return {k : len(v) for k, v in child_nodes.iteritems()}

    def mergers_per_level(self):
        """
        Levels at which mergers occur
        Returns {level : n_mergers}
        """
        merge_levels = [d[0] + 1 for d in self.mergers.values()]
        mergers = Counter(merge_levels)
        return mergers


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
        return fig, ax
        # fig.show()

