import matplotlib.pyplot as plt
from netCDF4 import Dataset
from homology.util.animations import make_animation
plt.rcParams['figure.figsize'] = [10, 7]

# This produces the animations for LES data, which will compare the 4 different datasets, always for the
# vertical velocity field. We generate a time sweep along an arbitrary vertical plane and a height sweep.

data_path = '/media/licon/1650BDF550BDDBA54/CBL3D_Data_LES/Shaofeng_Third/'
movies_path = '/home/licon/uni-koeln/tr32/notebooks/movies/vertical_wind_LES/'

fixed_time = 555
fixed_y = 25

################
# Dataset 1 ####
################

infile = data_path + 'SP1/wrfout_d01_2009-08-05_08-00-00'
dataset = Dataset(infile, 'r')
vardata = dataset.variables['W']

# Height sweep
plane_axes = (2, 3)
sweep_axis = 1
fixed_axis = 0
fixed_axis_value = fixed_time
outfile = movies_path + 'SP1_fix_t_sweep_z.gif'
tmp = make_animation(vardata, plane_axes, sweep_axis, fixed_axis, fixed_axis_value)
tmp.save(outfile, writer='imagemagick', fps=8)

# Time sweep
plane_axes = (1, 2)
sweep_axis = 0
fixed_axis = 3
fixed_axis_value = fixed_y
outfile = movies_path + 'SP1_fix_y_sweep_t.gif'
tmp = make_animation(vardata, plane_axes, sweep_axis, fixed_axis, fixed_axis_value)
tmp.save(outfile, writer='imagemagick', fps=4)
# close connection
dataset.close()

################
# Dataset 2 ####
################

infile = data_path + 'SP2/wrfout_d01_2009-08-05_08-00-00'
dataset = Dataset(infile, 'r')
vardata = dataset.variables['W']

# Height
plane_axes = (2, 3)
sweep_axis = 1
fixed_axis = 0
fixed_axis_value = fixed_time
outfile = movies_path + 'SP2_fix_t_sweep_z.gif'
tmp = make_animation(vardata, plane_axes, sweep_axis, fixed_axis, fixed_axis_value)
tmp.save(outfile, writer='imagemagick', fps=8)

# Time
plane_axes = (1, 2)
sweep_axis = 0
fixed_axis = 3
fixed_axis_value = fixed_y
outfile = movies_path + 'SP2_fix_y_sweep_t.gif'
tmp = make_animation(vardata, plane_axes, sweep_axis, fixed_axis, fixed_axis_value)
tmp.save(outfile, writer='imagemagick', fps=4)
# close connection
dataset.close()

################
# Dataset 3 ####
################

infile = data_path + 'SP3/wrfout_d01_2009-08-05_08-00-00'
dataset = Dataset(infile, 'r')
vardata = dataset.variables['W']

# Height
plane_axes = (2, 3)
sweep_axis = 1
fixed_axis = 0
fixed_axis_value = fixed_time
outfile = movies_path + 'SP3_fix_t_sweep_z.gif'
tmp = make_animation(vardata, plane_axes, sweep_axis, fixed_axis, fixed_axis_value)
tmp.save(outfile, writer='imagemagick', fps=8)

# Time
plane_axes = (1, 2)
sweep_axis = 0
fixed_axis = 3
fixed_axis_value = fixed_y
outfile = movies_path + 'SP3_fix_y_sweep_t.gif'
tmp = make_animation(vardata, plane_axes, sweep_axis, fixed_axis, fixed_axis_value)
tmp.save(outfile, writer='imagemagick', fps=4)
# close connection
dataset.close()

################
# Dataset 4 ####
################

infile = data_path + 'SP4/wrfout_d01_2009-08-05_08-00-00'
dataset = Dataset(infile, 'r')
vardata = dataset.variables['W']

# Height
plane_axes = (2, 3)
sweep_axis = 1
fixed_axis = 0
fixed_axis_value = fixed_time
outfile = movies_path + 'SP4_fix_t_sweep_z.gif'
tmp = make_animation(vardata, plane_axes, sweep_axis, fixed_axis, fixed_axis_value)
tmp.save(outfile, writer='imagemagick', fps=8)

# Time
plane_axes = (1, 2)
sweep_axis = 0
fixed_axis = 3
fixed_axis_value = fixed_y
outfile = movies_path + 'SP4_fix_y_sweep_t.gif'
tmp = make_animation(vardata, plane_axes, sweep_axis, fixed_axis, fixed_axis_value)
tmp.save(outfile, writer='imagemagick', fps=4)
# close connection
dataset.close()