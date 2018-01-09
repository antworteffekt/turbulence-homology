import numpy as np
from netCDF4 import Dataset
from homology.mergetree import MergeTree
import dill

# Open connection
dataset = Dataset('/media/licon/1650BDF550BDDBA54/CBL3D_Data_LES/Shaofeng_Third/SP1/wrfout_d01_2009-08-05_08-00-00', 'r')
# out_dir = '/media/licon/1650BDF550BDDBA54/mergetrees/'
out_dir = '/media/licon/1650BDF550BDDBA54/mergetrees/3d/'
# out_dir = '/home/licon/uni-koeln/tr32/tmp/'
# out_dir = '/home/licon/uni-koeln/tr32/tmp/3d/'
# Set coordinate value for the 2D plane
fixed_axis_value = 25

# t_range = np.arange(dataset.variables['W'].shape[0])
t_range = np.arange(500, 520)

for t in t_range:
    print "Timestep %d" % t
    tree = MergeTree()
    tree.build_merge_tree_3d(dataset.variables['W'][t, :, :, :], tol=0)
    outfile_name = out_dir + 'tree_%d.pickle' % t
    with open(outfile_name, 'w+') as f:
        dill.dump(tree, f)
    # print outfile_name

# Close connection
dataset.close()
