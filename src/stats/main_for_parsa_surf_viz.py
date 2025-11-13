from dfsio import readdfs, writedfs
import os
from surfproc import mean_curvature, view_patch_vtk, patch_color_attrib, smooth_patch
import numpy as np
from scipy.io import loadmat

# load mat file /home/ajoshi/Downloads/c.mat

cmat = loadmat('/home/ajoshi/Downloads/c.mat')

c=cmat['c']

lsurf = readdfs('/home/ajoshi/Projects/bfp/supp_data/bci32kleft.dfs')

c = c[:len(lsurf.vertices)]

l = type('', (), {})()  # create empty object

l.vertices = lsurf.vertices.copy()
l.faces = lsurf.faces.copy()
l.attributes = c.flatten()

lsurf.attributes = mean_curvature(lsurf)

# clamp curvature values for better visualization
# Apply sigmoid transformation to clamp values smoothly
lsurf.attributes = 1 / (1 + np.exp(-lsurf.attributes * 50))

lsurf=patch_color_attrib(lsurf, cmap='Greys',clim=[0,1.5])

l = patch_color_attrib(l, cmap='jet', clim=[np.min(c), np.max(c)])

lsurf.vColor[c.flatten()> 0.15,:] = l.vColor[c.flatten()>0.15,:]  # color vertices with curvature > 0.5 with c values 

# smooth surface 
lsurf = smooth_patch(lsurf, iterations=1000)


view_patch_vtk(lsurf, azimuth=100, elevation=180, roll=90,
               outfile='left1.png', show=1)

writedfs('left_curv.dfs', lsurf)


