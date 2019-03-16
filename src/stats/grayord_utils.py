""" This module contains helful utility functions for handling and visualizing gray ordinate file """
import os
import scipy as sp
import numpy as np
import scipy.io as spio
import sys
sys.path.append('src/stats')
from surfproc import view_patch_vtk, patch_color_attrib, smooth_surf_function, smooth_patch
from dfsio import readdfs, writedfs


def visdata_grayord(data, surf_name, out_dir, smooth_iter, colorbar_lim,
                    colormap):
    lsurf, rsurf = label_surf(data, colorbar_lim, smooth_iter, colormap)
    vis_save_surf(lsurf, rsurf, out_dir, surf_name)


def vis_grayord_sigcorr(pval, rval, surf_name, out_dir, smooth_iter):
    print('outputs will be written to directory: ' + out_dir)
    plsurf, prsurf = label_surf(pval, [0, 0.05], smooth_iter, 'jet_r')
    # If p value above .05 then make the surface grey
    plsurf.vColor[plsurf.attributes >= 0.05, :] = .5
    prsurf.vColor[prsurf.attributes >= 0.05, :] = .5
    vis_save_surf(plsurf, prsurf, out_dir, surf_name + 'pval_sig')
    print('output pvalues on surface')
    print('colorbar limits are 0 to 0.05; colorbar class is jet reverse')

    lsurf, rsurf = label_surf(rval, [-.5, .5], smooth_iter, 'jet')
    # If p value above .05 then make the surface grey
    lsurf.vColor[plsurf.attributes >= 0.05, :] = .5
    rsurf.vColor[prsurf.attributes >= 0.05, :] = .5
    vis_save_surf(lsurf, rsurf, out_dir, surf_name + 'rval_sig')
    print('output pvalues on surface')
    print('colorbar limits are -0.5 to +0.5; colorbar class is jet')


def vis_grayord_sigpval(pval, surf_name, out_dir, smooth_iter, bfp_path='.'):
    plsurf, prsurf = label_surf(
        pval, [0, 0.05], smooth_iter, 'jet_r', bfp_path=bfp_path)
    # If p value above .05 then make the surface grey
    plsurf.vColor[plsurf.attributes >= 0.05, :] = .5
    prsurf.vColor[prsurf.attributes >= 0.05, :] = .5
    vis_save_surf(
        plsurf, prsurf, out_dir, surf_name + 'pval_sig', bfp_path=bfp_path)


def label_surf(pval, colorbar_lim, smooth_iter, colormap, bfp_path='.'):
    lsurf = readdfs(os.path.join(bfp_path, 'supp_data/bci32kleft.dfs'))
    rsurf = readdfs(os.path.join(bfp_path, 'supp_data/bci32kright.dfs'))
    num_vert = lsurf.vertices.shape[0]
    lsurf.attributes = sp.zeros((lsurf.vertices.shape[0]))
    rsurf.attributes = sp.zeros((rsurf.vertices.shape[0]))
    #smooth surfaces
    lsurf = smooth_patch(lsurf, iterations=smooth_iter)

    rsurf = smooth_patch(rsurf, iterations=smooth_iter)

    # write on surface attributes
    lsurf.attributes = pval.squeeze()
    lsurf.attributes = lsurf.attributes[:num_vert]
    rsurf.attributes = pval.squeeze()
    rsurf.attributes = rsurf.attributes[num_vert:2 * num_vert]

    lsurf = patch_color_attrib(lsurf, clim=colorbar_lim, cmap=colormap)
    rsurf = patch_color_attrib(rsurf, clim=colorbar_lim, cmap=colormap)

    return lsurf, rsurf


def vis_save_surf(lsurf, rsurf, out_dir, surf_name, bfp_path='.'):
    # if label is zero, black out surface, attribute should be nan
    num_vert = lsurf.vertices.shape[0]
    lab = spio.loadmat(
        os.path.join(bfp_path, 'supp_data/USCBrain_grayordinate_labels.mat'))
    labs = lab['labels'].squeeze()
    labs = sp.float64(labs)
    lsurf.attributes[labs[:num_vert] == 0] = sp.nan
    rsurf.attributes[labs[num_vert:2 * num_vert] == 0] = sp.nan
    lsurf.vColor[sp.isnan(lsurf.attributes), :] = 0
    rsurf.vColor[sp.isnan(lsurf.attributes), :] = 0

    # Visualize left hemisphere
    view_patch_vtk(
        lsurf,
        azimuth=100,
        elevation=180,
        roll=90,
        outfile=out_dir + '/LeftLateral_' + surf_name + '.png',
        show=0)
    view_patch_vtk(
        lsurf,
        azimuth=-100,
        elevation=180,
        roll=-90,
        outfile=out_dir + '/LeftMedial_' + surf_name + '.png',
        show=0)
    # Visualize right hemisphere
    view_patch_vtk(
        rsurf,
        azimuth=-100,
        elevation=180,
        roll=-90,
        outfile=out_dir + '/RightLateral_' + surf_name + '.png',
        show=0)
    view_patch_vtk(
        rsurf,
        azimuth=100,
        elevation=180,
        roll=90,
        outfile=out_dir + '/RightMedial_' + surf_name + '.png',
        show=0)

    writedfs(out_dir + '/Right_' + surf_name + '.dfs', rsurf)
    writedfs(out_dir + '/Left_' + surf_name + '.dfs', lsurf)
