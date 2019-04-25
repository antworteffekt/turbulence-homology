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
                # print type(cell_parent)
                # print 'Cell: %s  --- parent: %s' % (it.multi_index, cell_parent)
                # Check for connectivity behind and below only - this works because
                # we assume a two-dimensional array structure
                # Periodicity should affect only the lookup of the cell "behind".
                if periodic:
                    behind = tuple(map(lambda x, y: x + y, it.multi_index, (0, -1)))
                    # Map the column index to its residue mod n_columns.
                    behind = (behind[0], behind[1] % X.shape[1])
                    if X[behind]:
                        self.uf.union(it.multi_index, behind, elder_rule=False)
                else:
                    behind = tuple(map(lambda x, y: x + y, it.multi_index, (0, -1)))
                    if not any(n < 0 for n in behind) and X[behind]:
                        self.uf.union(it.multi_index, behind, elder_rule=False)

                below = tuple(map(lambda x, y: x + y, it.multi_index, (-1, 0)))
                if not any(n < 0 for n in below) and X[below]:
                    # Don't do the union here, just keep track of the involved
                    # elements instead.
                    # Get cell's parent in this level
                    cell_parent = self.uf.find(it.multi_index)
                    # Get root of cell below
                    below_parent = self.uf.find(below)
                    if cell_parent not in below_joins:
                        below_joins[cell_parent] = [tuple(below_parent)]
                    else:
                        below_joins[cell_parent].append(below_parent)

            it.iternext()
            if it.finished or it.multi_index[0] != current_level:
                for level_parent, below_parent in below_joins.items():
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
    # @profile
    def build_merge_tree_3d(self, X, periodic=False, return_graph=False, iter_dim=0, flip=False,
                            clear_uf=False, **kwargs):
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
                    behind = tuple(map(lambda x, y: x + y, it.multi_index, (0, 0, -1)))
                    behind = (behind[0], behind[1], behind[2] % X.shape[2])
                    if X[behind]:
                        self.uf.union(it.multi_index, behind, elder_rule=False)

                    side = tuple(map(lambda a, b: a + b, it.multi_index, (0, -1, 0)))
                    side = (side[0], side[1] % X.shape[1], side[2])
                    if X[side]:
                        self.uf.union(it.multi_index, side, elder_rule=False)
                else:
                    # Check for connectivity behind, below, and to one side - this assumes 
                    # six connectivity
                    behind = tuple(map(lambda x, y: x + y, it.multi_index, (0, 0, -1)))
                    if not any(n < 0 for n in behind) and X[behind]:
                        self.uf.union(it.multi_index, behind, elder_rule=False)

                    side = tuple(map(lambda a, b: a + b, it.multi_index, (0, -1, 0)))
                    if not any(n < 0 for n in side) and X[side]:
                        self.uf.union(it.multi_index, side, elder_rule=False)

                below = tuple(map(lambda x, y: x + y, it.multi_index, (-1, 0, 0)))
                if not any(n < 0 for n in below) and X[below]:
                    # Don't do the union here, just keep track of the involved elements instead.
                    # Get cell's parent in this level
                    cell_parent = self.uf.find(it.multi_index)
                    # Get root of cell below
                    below_parent = self.uf.find(below)
                    if cell_parent not in below_joins:
                        below_joins[cell_parent] = [tuple(below_parent)]
                    else:
                        below_joins[cell_parent].append(below_parent)

            it.iternext()
            if it.finished or it.multi_index[0] != current_level:
                # Level empty, store merger information
                for level_parent, below_parent in below_joins.items():
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

                # print current_level
                # print self.uf.objects_to_num
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
                            # Update values to new parent
                            children = [k for k, v in self.uf.parent_pointers.items() if v == self.uf.objects_to_num[cell]]
                            # print 'updating %d values ' % len(children)
                            for k in children:
                                self.uf.parent_pointers[k] = self.uf.objects_to_num[new_parent]
                            # if cell == (3, 199, 85):
                            #     print 'YOOT'
                            #     print 'cell - new parent', cell, new_parent
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
                # Delete superfluous information
                if clear_uf and current_level >= 2:
                    # print "*************\nLevel %d finished" % current_level
                    # elem = (3, 199, 85)
                    # try:
                    z_2 = [self.uf.num_to_objects[v] for k, v in self.uf.objects_to_num.items() if k[0] == current_level-2]
                    # if (3, 199, 85) in z_2:
                    #     print "WTF"
                    # except KeyError:
                    #     print current_level
                    #     print z_2
                    #     print self.uf.num_to_objects
                    #     raise
                    # print "Level z-2 has %d elements" % len(z_2)
                    # print "finding non roots, level %d" % current_level
                    # print set([k[0] for k in z_2])
                    # if current_level >= 3:
                        # print 'num. children of (3, 199, 85): ', len([k for k, v in self.uf.parent_pointers.items() if v == self.uf.objects_to_num[(3, 199, 85)]])
                        # print 'identifier: ', self.uf.objects_to_num[elem], 'element: ', elem
                        # print ' *** parent: ***', self.uf.find(elem)
                    non_roots = [x for x in z_2 if x != self.uf.find(x)]
                    # roots = [x for x in z_2 if x == self.uf.find(x)]
                    # print '%d roots found ' % len(roots), roots
                    # print "Non-root elements: %d" % len(non_roots)
                    # print non_roots
                    for t in non_roots:
                        # self.uf.find(t)
                        identifier = self.uf.objects_to_num.pop(t)
                        # if identifier == 568144:
                        #     print 'element deleted at level ', current_level
                        #     print identifier, ' --- ', t
                        self.uf.num_to_objects.pop(identifier)
                        self.uf.parent_pointers.pop(identifier)
                    # print '************'
        # Clear the last two levels
        if clear_uf:
            levels = X.shape[0]-1, X.shape[0]-2
            for z in levels:
                zz = [self.uf.num_to_objects[v] for k, v in self.uf.objects_to_num.items() if k[0] == z]
                non_roots = [x for x in zz if x != self.uf.find(x)]
                for t in non_roots:
                    identifier = self.uf.objects_to_num.pop(t)
                    self.uf.num_to_objects.pop(identifier)
                    self.uf.parent_pointers.pop(identifier)
        if return_graph:
            return G

    def structure_size(self):
        """
        Count the number of cells in each connected structure.
        Returns {parent : #(cells)}
        """
        # First invert the dict parent_pointers
        child_nodes = {}
        for k, v in self.uf.parent_pointers.items():
            keys = child_nodes.setdefault(v, [])
            keys.append(k)
        return {k : len(v) for k, v in child_nodes.items()}

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
        for k, d in self.components_level.items():
            ax.vlines(k[1], d[0], d[-1])
        for k, d in self.mergers.items():
            xmin = min(k[1], d[1][1])
            xmax = max(k[1], d[1][1])
            ax.hlines(d[0] + 1, xmin, xmax)
            ax.vlines(k[1], d[0], d[0] + 1)
            ax.plot(d[1][1], d[0] + 1, 'r.')
        # ax.invert_xaxis()
        return fig, ax
        # fig.show()


def find_root(n, tree, root_list):
    if n not in root_list:
        return find_root(tree.uf.parent_pointers[n], tree, root_list)
    else:
        return n


def update_parent_pointers(tree):

    A = set(tree.uf.parent_pointers.values())
    B = set(tree.uf.num_weights.keys())

    diff = A - B
    if len(diff) == 0:
        return
    for n in diff:
        map(tree.uf.find,
            [tree.uf.num_to_objects[k] for k, v in tree.uf.parent_pointers.items()
             if v == n])


def calculate_wsquared_level(X, tree):

    # Number of height levels
    n_levels = X.shape[0]
    levels = range(n_levels)
    components_values = {k: [0]*n_levels for k in set(tree.uf.num_weights.keys())}

    for l in levels:
        level_cells = {k: tree.uf.parent_pointers[v]
                       for k, v in tree.uf.objects_to_num.items() if k[0] == l}
        for k, v in level_cells.items():
            try:
                components_values[v][l] += X[k] ** 2
            except KeyError:
                key_cell = tree.uf.num_to_objects[v]
                parent = tree.uf.objects_to_num[tree.uf.find(key_cell)]
                components_values[parent][l] += X[k] ** 2

    return components_values


def calculate_wb_level(W, T, tree):

    # Number of height levels
    n_levels = W.shape[0]
    levels = range(n_levels)
    components_values = {k: [0]*n_levels for k in set(tree.uf.num_weights.keys())}

    for l in levels:
        level_cells = {k: tree.uf.parent_pointers[v]
                       for k, v in tree.uf.objects_to_num.items() if k[0] == l}
        for k, v in level_cells.items():
            try:
                components_values[v][l] += W[k] * T[k]
            except KeyError:
                key_cell = tree.uf.num_to_objects[v]
                parent = tree.uf.objects_to_num[tree.uf.find(key_cell)]
                components_values[parent][l] += W[k] * T[k]

    return components_values


def calculate_variables_level(W, T, tree):

    # Number of height levels:
    n_levels = W.shape[0]
    levels = range(n_levels)
    # w_squared_level = {k: np.zeros(n_levels) for k in set(tree.uf.num_weights.keys())}
    # wt_level = {k: np.zeros(n_levels) for k in set(tree.uf.num_weights.keys())}
    # sizes_level = {k: np.zeros(n_levels) for k in set(tree.uf.num_weights.keys())}
    w_squared_level = {k: np.zeros(n_levels) for k in set(tree.uf.parent_pointers.values())}
    wt_level = {k: np.zeros(n_levels) for k in set(tree.uf.parent_pointers.values())}
    sizes_level = {k: np.zeros(n_levels) for k in set(tree.uf.parent_pointers.values())}

    for l in levels:
        level_cells = {k: tree.uf.parent_pointers[v]
                       for k, v in tree.uf.objects_to_num.items() if k[0] == l}
        counts = Counter(level_cells.values())
        # for k, v in level_cells.items():
        #     try:
        #         w_squared_level[v][l]  += W[k] ** 2
        #         wt_level[v][l] += W[k] * T[k]
        #         sizes_level[v][l] = counts[v]
        #     except KeyError:
        #         key_cell = tree.uf.num_to_objects[v]
        #         parent = tree.uf.objects_to_num[tree.uf.find(key_cell)]
        #         w_squared_level[parent][l] += W[k] ** 2
        #         wt_level[parent][l] += W[k] * T[k]
        #         sizes_level[parent][l] = counts[v]
        for k, v in counts.items():
            sizes_level[k][l] = v
        for k, v in level_cells.items():
            w_squared_level[v][l]  += W[k] ** 2
            wt_level[v][l] += W[k] * T[k]

    A = set(tree.uf.parent_pointers.values())
    B = set(tree.uf.num_weights.keys())
    diff = A - B
    update_map = {k: find_root(k, tree, B) for k in diff}
    w_squared = {k: np.zeros(n_levels) for k in tree.uf.num_weights.keys()}
    wt = {k: np.zeros(n_levels) for k in tree.uf.num_weights.keys()}
    sizes = {k: np.zeros(n_levels) for k in tree.uf.num_weights.keys()}
    
    # for k, v in sizes_level.items():
    #     if k in sizes.keys():
    #         sizes[k] += v
    #     else:
    #         sizes[update_map[k]] += v
    for k, v in sizes_level.items():
        if k in sizes.keys():
            sizes[k] += v
        else:
            sizes[update_map[k]] += v
    for k, v in w_squared_level.items():
        if k in w_squared.keys():
            w_squared[k] += v
        else:
            w_squared[update_map[k]] += v
    for k, v in wt_level.items():
        if k in wt.keys():
            wt[k] += v
        else:
            wt[update_map[k]] += v

    return (w_squared, wt, sizes)


def calculate_sizes_level(X, tree):

    # Number of height levels
    n_levels = X.shape[0]
    levels = range(n_levels)
    sizes = {k: np.zeros(n_levels) for k in set(tree.uf.parent_pointers.values())}

    for l in levels:
        level_cells = {k: tree.uf.parent_pointers[v]
                       for k, v in tree.uf.objects_to_num.items() if k[0] == l}
        counts = Counter(level_cells.values())
        for k, v in counts.items():
            sizes[k][l] = v

    A = set(tree.uf.parent_pointers.values())
    B = set(tree.uf.num_weights.keys())
    diff = A - B
    update_map = {k: find_root(k, tree, B) for k in diff}
    components_sizes = {k: np.zeros(n_levels) for k in tree.uf.num_weights.keys()}
    for k, v in sizes.items():
        if k in components_sizes.keys():
            components_sizes[k] += v
        else:
            components_sizes[update_map[k]] += v

    return components_sizes
