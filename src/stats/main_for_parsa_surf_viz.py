from dfsio import readdfs, writedfs
import os
from surfproc import mean_curvature, view_patch_vtk, patch_color_attrib, smooth_patch
import numpy as np
from scipy.io import loadmat

# load mat file /home/ajoshi/Downloads/c.mat

cmat = loadmat('/home/ajoshi/Downloads/c.mat')

lsurf = readdfs('/home/ajoshi/Projects/bfp/supp_data/bci32kleft.dfs')

lsurf.attributes = mean_curvature(lsurf)

# clamp curvature values for better visualization
# Apply sigmoid transformation to clamp values smoothly
lsurf.attributes = 1 / (1 + np.exp(-lsurf.attributes * 50))

lsurf=patch_color_attrib(lsurf, cmap='Greys',clim=[0,1.5])


# smooth surface 
lsurf = smooth_patch(lsurf, iterations=1000)


view_patch_vtk(lsurf, azimuth=100, elevation=180, roll=90,
               outfile='left1.png', show=1)

writedfs('left_curv.dfs', lsurf)


