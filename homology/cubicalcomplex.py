import numpy as np
import matplotlib.pyplot as plt


class CubicalComplex:

    def __init__(self, dim=2):
        self.cubes = {}
        self.anchor_points = None
        if dim > 2:
            raise NotImplementedError("Dimensions > 2 not yet implemented!")
        self.dim = dim
        for d in range(dim + 1):
            self.cubes[d] = None

    def from_anchor_points(self, points):

        self.anchor_points = points
        if points.shape[1] != self.dim:
            raise ValueError("Dimensions don't agree.")
        # Add the 0-cells
        self.cubes[0] = points
        corners = []
        for p in points:
            corners.extend(
                [(p[0] + 1, p[1]), (p[0], p[1] + 1), (p[0] + 1, p[1] + 1)])
        self.cubes[0] = np.unique(np.concatenate(
            [self.cubes[0], np.array(corners)]), axis=0)
        # Add the 1-cells spawning from each anchor point
        self.cubes[1] = self._spawn_edges()
        # 2-cells are not added explicitly at this point

    def _spawn_edges(self):
        edges = []
        for p in self.anchor_points:
            edges.extend([((p[0], p[1]), (p[0], p[1] + 1)),
                          ((p[0], p[1]), (p[0] + 1, p[1])),
                          ((p[0], p[1] + 1), (p[0] + 1, p[1] + 1)),
                          ((p[0] + 1, p[1]), (p[0] + 1, p[1] + 1))])
        return list(set(edges))

    def plot(self, ax=None, alpha=0.3, **kwargs):

        if ax is None:
            f = plt.figure()
            sp = f.add_subplot(111, aspect='equal')
        else:
            sp = ax
        sp.axis('off')
        sp.plot(self.cubes[0][:, 0], self.cubes[0][
                :, 1], 'o', color='black', **kwargs)
        for l in self.cubes[1]:
            sp.plot((l[0][0], l[1][0]), (l[0][1], l[1][1]),
                    color='black', lw=1.5)
        for p in self.anchor_points:
            rect = plt.Rectangle((p[0], p[1]), 1, 1, color='black', alpha=alpha)
            sp.add_patch(rect)
