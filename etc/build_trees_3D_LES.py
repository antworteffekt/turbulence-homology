import numpy as np
from netCDF4 import Dataset
from homology.mergetree import MergeTree
import dill

# Dataset SP1
# Open connection
dataset = Dataset('/media/licon/1650BDF550BDDBA54/CBL3D_Data_LES/Shaofeng_Third/SP1/wrfout_d01_2009-08-05_08-00-00', 'r')
# out_dir = '/media/licon/1650BDF550BDDBA54/mergetrees/'
out_dir = '/media/licon/1650BDF550BDDBA54/mergetrees/3d/vertical_velocity/SP1/'
# out_dir = '/home/licon/uni-koeln/tr32/tmp/3d/'

# Threshold value to use in construction of binary array
threshold = 0.01

# Timestep range to process
t_range = np.arange(400, 721)

for t in t_range:
    print "Timestep %d" % t
    tree = MergeTree()
    tree.build_merge_tree_3d(dataset.variables['W'][t, :, :, :], tol=threshold)
    outfile_name = out_dir + 'tree_%d.pickle' % t
    with open(outfile_name, 'w+') as f:
        dill.dump(tree, f)

# Close connection
dataset.close()
