import numpy as np
import matplotlib.pyplot as plt

from matplotlib import animation, cm
from homology.util.visualize import shiftedColorMap

# from IPython.display import HTML

# from netCDF4 import Dataset

# plt.rcParams['figure.figsize'] = [10, 7]


def make_animation(X, plane_axes, sweep_axis, fixed_axis, fixed_axis_value):
    """
    Returns animation object from specified data
    X                 A netCDF4 Dataset object
    plane_axes
    sweep_axis
    fixed_axis
    fixed_axis_value
    """

    # Create the array indexes
    obj = [0, 0, 0, 0]
    i = 0
    for idx in plane_axes:
        obj[idx] = slice(None, None, None)
    obj[fixed_axis] = fixed_axis_value
    obj[sweep_axis] = i

    def get_middle(xx):
        return 1 - (float(np.amax(xx)) / (np.amax(xx) + abs(np.amin(xx))))

    def init():
        global fig, ax, im, tx
        fig = plt.figure()
        ax = plt.axes()
        idx = list(obj)
        idx[sweep_axis] = slice(None, None, None)
        middle = get_middle(X[idx])
        # print obj
        im = ax.imshow(X[obj], cmap=shiftedColorMap(cm.seismic, midpoint=middle),
                       interpolation='none', aspect='auto',
                       vmin=np.amin(X[idx]), vmax=np.amax(X[idx]))
        ax.invert_yaxis()
        cb = fig.colorbar(im)
        tx = ax.set_title('%s = %d' % (X.dimensions[sweep_axis], i))
        return

    def animate(n):
        # update indexes
        obj[sweep_axis] = n
        # vmax = np.max(X[obj])
        # vmin = np.min(X[obj])
        im.set_data(X[obj])
        # im.set_clim(vmin, vmax)
        tx.set_text('%s = %d' % (X.dimensions[sweep_axis], n))
        return

    init()
    anim = animation.FuncAnimation(fig, animate, frames=X.shape[sweep_axis], interval=100, blit=False)
    return anim

    # TEST PLOT
    # fig, ax = plt.subplots()
    # # interp = kwargs.get('interpolation', 'none')
    # # colors = kwargs.get('colormap', 'seismic')
    # img0 = ax.imshow(X[obj], cmap='Blues', interpolation='none')
    # fig.colorbar(img0, ax=ax, fraction=0.022, pad=0.01)
    # ax.invert_yaxis()
    # tx = ax.set_title('%s = %d' % (X.dimensions[sweep_axis], i))
    # fig.tight_layout()
    # fig.show()

def make_animation_subset_levels(X, fixed_axes, fixed_value_1, fixed_value_2,
                                 filtration_size):
    """
    Returns animation object from specified data
    X                 A netCDF4 Dataset object
    fixed_axes
    fixed_value_1
    fixed_value_2
    filtration_size
    """
    # Create the array indexes
    obj = [slice(None, None, None)] * 4
    obj[fixed_axes[0]] = fixed_value_1
    obj[fixed_axes[1]] = fixed_value_2
    # print obj

    # Create sequence of threshold values
    thresholds = np.linspace(start=np.amin(X[obj]), stop=np.amax(X[obj]), num=filtration_size)
    # print thresholds
    # TEST PLOT
    # fig, ax = plt.subplots()
    # # interp = kwargs.get('interpolation', 'none')
    # # colors = kwargs.get('colormap', 'seismic')
    # img0 = ax.imshow(X[obj], cmap='Blues', interpolation='none')
    # fig.colorbar(img0, ax=ax, fraction=0.022, pad=0.01)
    # ax.invert_yaxis()
    # # tx = ax.set_title('%s = %d' % (X.dimensions[sweep_axis], i))
    # fig.tight_layout()
    # fig.show()

    # def get_middle(xx):
    #     return 1 - (float(np.amax(xx)) / (np.amax(xx) + abs(np.amin(xx))))

    def init():
        global fig, ax, im, tx
        fig = plt.figure()
        ax = plt.axes()
        # idx = list(obj)
        # idx[sweep_axis] = slice(None, None, None)
        # middle = get_middle(X[idx])
        # print obj
        im = ax.imshow(X[obj] < thresholds[2], cmap='Blues',#cmap=shiftedColorMap(cm.seismic, midpoint=middle),
                       interpolation='none', aspect='auto')
                       # vmin=np.amin(X[idx]), vmax=np.amax(X[idx]))
        ax.invert_yaxis()
        # cb = fig.colorbar(im)
        # tx = ax.set_title('%s = %d' % (X.dimensions[sweep_axis], i))
        return

    def animate(n):
        # update indexes
        # obj[sweep_axis] = n
        # vmax = np.max(X[obj])
        # vmin = np.min(X[obj])
        im.set_data(X[obj] < thresholds[n])
        # im.set_clim(vmin, vmax)
        # tx.set_text('%s = %d' % (X.dimensions[sweep_axis], n))
        return

    init()
    anim = animation.FuncAnimation(fig, animate, frames=np.arange(filtration_size), interval=100, blit=False)
    return anim
