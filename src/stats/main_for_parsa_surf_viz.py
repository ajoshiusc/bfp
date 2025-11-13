from dfsio import readdfs, writedfs
import os
from surfproc import mean_curvature, view_patch_vtk, patch_color_attrib, smooth_patch
import numpy as np
from scipy.io import loadmat

# load mat file /home/ajoshi/Downloads/c.mat

cmat = loadmat("/home/ajoshi/Downloads/c.mat")

c = cmat["c"]

lsurf = readdfs("/home/ajoshi/Projects/bfp/supp_data/bci32kleft.dfs")

c = c[: len(lsurf.vertices)]

l = type("", (), {})()  # create empty object

l.vertices = lsurf.vertices.copy()
l.faces = lsurf.faces.copy()
l.attributes = c.flatten()

lsurf.attributes = mean_curvature(lsurf)

# clamp curvature values for better visualization
# Apply sigmoid transformation to clamp values smoothly
lsurf.attributes = 1 / (1 + np.exp(-lsurf.attributes * 50))

lsurf = patch_color_attrib(lsurf, cmap="Greys", clim=[0, 1.5])
# smooth surface
lsurf = smooth_patch(lsurf, iterations=1000)

writedfs("left_curv_orig.dfs", lsurf)

l = patch_color_attrib(l, cmap="jet", clim=[np.min(c), np.max(c)])

lsurf.vColor[np.abs(c.flatten()) > 0.15, :] = l.vColor[
    np.abs(c.flatten()) > 0.15, :
]  # color vertices with curvature > 0.5 with c values


view_patch_vtk(lsurf, azimuth=100, elevation=180, roll=90, outfile="left1.png", show=1)

writedfs("left_curv.dfs", lsurf)

# add the code for right hemisphere as well
rsurf = readdfs("/home/ajoshi/Projects/bfp/supp_data/bci32kright.dfs")
rsurf.attributes = mean_curvature(rsurf)
rsurf.attributes = 1 / (1 + np.exp(-rsurf.attributes * 50))
rsurf = patch_color_attrib(rsurf, cmap="Greys", clim=[0, 1.5])
# smooth surface
rsurf = smooth_patch(rsurf, iterations=1000)
writedfs("right_curv_orig.dfs", rsurf)

c = cmat["c"]
c = c[len(lsurf.vertices) : len(lsurf.vertices) + len(rsurf.vertices)]
r = type("", (), {})()  # create empty object
r.vertices = rsurf.vertices.copy()
r.faces = rsurf.faces.copy()
r.attributes = c.flatten()
r = patch_color_attrib(r, cmap="jet", clim=[np.min(c), np.max(c)])
rsurf.vColor[np.abs(c.flatten()) > 0.15, :] = r.vColor[np.abs(c.flatten()) > 0.15, :]  # color vertices with curvature > 0.5 with c values
view_patch_vtk(rsurf, azimuth=100, elevation=180, roll=90, outfile="right1.png", show=1)
writedfs("right_curv.dfs", rsurf)
