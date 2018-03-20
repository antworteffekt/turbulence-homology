import numpy as np
from netCDF4 import Dataset
from homology.mergetree import MergeTree
import dill

# Dataset SP1
# Open connection
dataset = Dataset('/media/licon/1650BDF550BDDBA54/CBL3D_Data_LES/Shaofeng_Third/SP1/wrfout_d01_2009-08-05_08-00-00', 'r')
out_dir = '/home/licon/uni-koeln/tr32/tmp/SP1/'
# Set coordinate value for the 2D plane
fixed_axis_value = 25
t_range = np.arange(dataset.variables['W'].shape[0])

for t in t_range:
    print "Timestep %d" % t
    tree = MergeTree()
    tree.build_merge_tree(dataset.variables['W'][t, :, :, fixed_axis_value], tol=0.01)
    outfile_name = out_dir + 'tree_%d.pickle' % t
    with open(outfile_name, 'w+') as f:
        dill.dump(tree, f)

# Close connection
dataset.close()

##############################################################################################################

# Dataset SP2
# Open connection
dataset = Dataset('/media/licon/1650BDF550BDDBA54/CBL3D_Data_LES/Shaofeng_Third/SP2/wrfout_d01_2009-08-05_08-00-00', 'r')
out_dir = '/home/licon/uni-koeln/tr32/tmp/SP2/'
# Set coordinate value for the 2D plane
fixed_axis_value = 25
t_range = np.arange(dataset.variables['W'].shape[0])

for t in t_range:
    print "Timestep %d" % t
    tree = MergeTree()
    tree.build_merge_tree(dataset.variables['W'][t, :, :, fixed_axis_value], tol=0.01)
    outfile_name = out_dir + 'tree_%d.pickle' % t
    with open(outfile_name, 'w+') as f:
        dill.dump(tree, f)

# Close connection
dataset.close()

##############################################################################################################

# Dataset SP3
# Open connection
dataset = Dataset('/media/licon/1650BDF550BDDBA54/CBL3D_Data_LES/Shaofeng_Third/SP3/wrfout_d01_2009-08-05_08-00-00', 'r')
out_dir = '/home/licon/uni-koeln/tr32/tmp/SP3/'
# Set coordinate value for the 2D plane
fixed_axis_value = 25
t_range = np.arange(dataset.variables['W'].shape[0])

for t in t_range:
    print "Timestep %d" % t
    tree = MergeTree()
    tree.build_merge_tree(dataset.variables['W'][t, :, :, fixed_axis_value], tol=0.01)
    outfile_name = out_dir + 'tree_%d.pickle' % t
    with open(outfile_name, 'w+') as f:
        dill.dump(tree, f)

# Close connection
dataset.close()

##############################################################################################################

# Dataset SP4
# Open connection
dataset = Dataset('/media/licon/1650BDF550BDDBA54/CBL3D_Data_LES/Shaofeng_Third/SP4/wrfout_d01_2009-08-05_08-00-00', 'r')
out_dir = '/home/licon/uni-koeln/tr32/tmp/SP4/'
# Set coordinate value for the 2D plane
fixed_axis_value = 25
t_range = np.arange(dataset.variables['W'].shape[0])

for t in t_range:
    print "Timestep %d" % t
    tree = MergeTree()
    tree.build_merge_tree(dataset.variables['W'][t, :, :, fixed_axis_value], tol=0.01)
    outfile_name = out_dir + 'tree_%d.pickle' % t
    with open(outfile_name, 'w+') as f:
        dill.dump(tree, f)

# Close connection
dataset.close()
