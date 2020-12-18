""" This module contains helful utility functions for handling and visualizing gray ordinate file """
import sys
import scipy.io as spio
import numpy as np
import scipy as sp
import os
from scipy.io import loadmat
from nilearn.image import load_img, new_img_like
from os.path import join
from dfsio import readdfs, writedfs
from surfproc import patch_color_attrib, smooth_surf_function
from scipy.ndimage import gaussian_filter

try:
    import vtk
    VTK_INSTALLED = 1
except ImportError as e:
    VTK_INSTALLED = 0

sys.path.append('src/stats')

if VTK_INSTALLED:
    from surfproc import view_patch_vtk, smooth_patch


FSL_PATH = '/usr/share/fsl/5.0'


def visdata_grayord(data,
                    surf_name,
                    out_dir,
                    smooth_iter,
                    colorbar_lim,
                    colormap,
                    save_png,
                    bfp_path='.',
                    fsl_path=FSL_PATH):
    lsurf, rsurf = label_surf(data,
                              colorbar_lim,
                              smooth_iter,
                              colormap,
                              bfp_path=bfp_path)
    save2surfgord(lsurf, rsurf, out_dir, surf_name, bfp_path, save_png)
    save2volgord(data, out_dir, surf_name, bfp_path)


def vis_grayord_sigcorr(pval, rval, sig_alpha, surf_name, out_dir, smooth_iter,
                        save_png, bfp_path, fsl_path):

    print('outputs will be written to directory: ' + out_dir)
    save2volgord(pval, out_dir, surf_name + '_pval', bfp_path, fsl_path)
    save2volgord(rval, out_dir, surf_name + '_rval', bfp_path, fsl_path)

    plsurf, prsurf = label_surf(pval, [0, sig_alpha], smooth_iter, 'jet_r',
                                bfp_path)
    # If p value above .05 then make the surface grey
    plsurf.vColor[plsurf.attributes > sig_alpha, :] = .5
    prsurf.vColor[prsurf.attributes > sig_alpha, :] = .5
    save2surfgord(plsurf, prsurf, out_dir, surf_name + '_pval_sig', bfp_path,
                  bool(save_png))
    print('output pvalues on surface')
    print('colorbar limits are 0 to ' + str(sig_alpha) +
          '; colorbar class is jet reverse')

    lsurf, rsurf = label_surf(rval, [-.5, .5], smooth_iter, 'jet', bfp_path)
    # If p value above .05 then make the surface grey
    lsurf.vColor[plsurf.attributes > sig_alpha, :] = .5
    rsurf.vColor[prsurf.attributes > sig_alpha, :] = .5
    save2surfgord(lsurf, rsurf, out_dir, surf_name + '_rval_sig', bfp_path,
                  bool(save_png))
    print('output rvalues on surface')
    print('colorbar limits are -0.5 to +0.5; colorbar class is jet')


def vis_grayord_sigpval(pval,
                        sig_alpha,
                        surf_name,
                        out_dir,
                        smooth_iter,
                        bfp_path,
                        fsl_path=FSL_PATH,
                        save_png=True):
    '''    save2volgord(pval,
                    out_dir,
                    surf_name + '_pval_sig',
                    bfp_path,
                    fsl_path=fsl_path,
                    default_value=0.5)
    '''

    save2volgord_bci((pval < sig_alpha)*(sig_alpha - pval)/sig_alpha,
                     out_dir,
                     surf_name + '_pval_sig',
                     bfp_path,
                     fsl_path=fsl_path)

    plsurf, prsurf = label_surf(pval, [0, sig_alpha],
                                smooth_iter,
                                'jet_r',
                                bfp_path=bfp_path)
    # If p value above .05 then make the surface grey
    plsurf.vColor[plsurf.attributes >= sig_alpha, :] = .5
    prsurf.vColor[prsurf.attributes >= sig_alpha, :] = .5
    save2surfgord(plsurf,
                  prsurf,
                  out_dir,
                  surf_name + 'pval_sig',
                  bfp_path=bfp_path,
                  save_png=save_png)


def label_surf(pval, colorbar_lim, smooth_iter, colormap, bfp_path='.'):
    lsurf = readdfs(os.path.join(bfp_path, 'supp_data/bci32kleft.dfs'))
    rsurf = readdfs(os.path.join(bfp_path, 'supp_data/bci32kright.dfs'))
    num_vert = lsurf.vertices.shape[0]
    lsurf.attributes = sp.zeros((lsurf.vertices.shape[0]))
    rsurf.attributes = sp.zeros((rsurf.vertices.shape[0]))

    if VTK_INSTALLED:
        # smooth surfaces
        lsurf = smooth_patch(lsurf, iterations=int(smooth_iter))
        rsurf = smooth_patch(rsurf, iterations=int(smooth_iter))
    else:
        print('VTK is not installed, surface will not be smoothed')

    # write on surface attributes
    lsurf.attributes = pval[:num_vert]  # .squeeze()
    #   lsurf.attributes = lsurf.attributes[:num_vert]
    rsurf.attributes = pval[num_vert:2 * num_vert]  # .squeeze()
    #    rsurf.attributes = rsurf.attributes[num_vert:2 * num_vert]

    lsurf = patch_color_attrib(lsurf, clim=colorbar_lim, cmap=colormap)
    rsurf = patch_color_attrib(rsurf, clim=colorbar_lim, cmap=colormap)

    return lsurf, rsurf


def save2volbord_bci(data, outfile, bfp_path='.', smooth_std=0):
    '''Save output to brainordinates'''
    a = loadmat(join(bfp_path, 'supp_data', 'bord_ind.mat'))
    v = load_img(join(bfp_path, 'supp_data',
                      'BCI-DNI_brain.pvc.frac.3mm.nii.gz'))

    img_dat = np.zeros(v.shape)
    img_dat[np.unravel_index(a['ind'], img_dat.shape,
                             order='F')] = data[:, None]

    if smooth_std > 0:
        img_dat = gaussian_filter(img_dat, sigma=(
            smooth_std, smooth_std, smooth_std), order=0)

    v2 = new_img_like(v, img_dat)

    v2.to_filename(outfile)


def save2volgord_bci(data, out_dir, vol_name, bfp_path='.', fsl_path=FSL_PATH, default_value=0, bst_path='/home/ajoshi/BrainSuite19b'):

    vol = load_img(
        join(bst_path, 'svreg', 'BCI-DNI_brain_atlas', 'BCI-DNI_brain.nii.gz'))
    a = loadmat(join(bfp_path, 'supp_data', 'bci_grayordinates_vol_ind.mat'))

    ind = a['bci_vol_ind'] - 1
    gordvol = np.zeros(vol.shape) + default_value

    val_gind = ~np.isnan(ind)
    ind2 = np.int32(ind[val_gind])
    ind2 = np.unravel_index(np.int32(ind2), vol.shape, order='F')
    gordvol[ind2] = data[val_gind.squeeze()]

    grod = new_img_like(vol, gordvol)
    grod.set_data_dtype(np.float64)
    outfile = join(out_dir, vol_name + '.nii.gz')
    grod.to_filename(outfile)


def save2volgord(data, out_dir, vol_name, bfp_path='.', fsl_path=FSL_PATH, default_value=0):

    mni2mm = load_img(join(fsl_path, 'data/standard', 'MNI152_T1_2mm.nii.gz'))
    a = loadmat(join(bfp_path, 'supp_data', 'MNI2mm_gord_vol_coord.mat'))

    ind = ~np.isnan(a['voxc']).any(axis=1)
    voxc = np.int16(a['voxc'] -
                    1)  # subtract 1 to convert from MATLAB to Python indexing
    gordvol = np.zeros(mni2mm.shape) + default_value

    gordvol[voxc[ind, 0], voxc[ind, 1], voxc[ind, 2]] = data[ind]
    grod = new_img_like(mni2mm, gordvol)
    grod.set_data_dtype(np.float64)
    outfile = join(out_dir, vol_name + '_MNI2mm.nii.gz')
    grod.to_filename(outfile)


def save2surfgord(lsurf,
                  rsurf,
                  out_dir,
                  surf_name,
                  bfp_path='.',
                  save_png=True):
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

    writedfs(out_dir + '/Right_' + surf_name + '.dfs', rsurf)
    writedfs(out_dir + '/Left_' + surf_name + '.dfs', lsurf)

    if VTK_INSTALLED == 0:
        print('VTK is not installed, screenshots will not be saved.')
        save_png = False

    if save_png == True:
        # Visualize left hemisphere
        view_patch_vtk(lsurf,
                       azimuth=100,
                       elevation=180,
                       roll=90,
                       outfile=out_dir + '/LeftLateral_' + surf_name + '.png',
                       show=0)
        view_patch_vtk(lsurf,
                       azimuth=-100,
                       elevation=180,
                       roll=-90,
                       outfile=out_dir + '/LeftMedial_' + surf_name + '.png',
                       show=0)
        # Visualize right hemisphere
        view_patch_vtk(rsurf,
                       azimuth=-100,
                       elevation=180,
                       roll=-90,
                       outfile=out_dir + '/RightLateral_' + surf_name + '.png',
                       show=0)
        view_patch_vtk(rsurf,
                       azimuth=100,
                       elevation=180,
                       roll=90,
                       outfile=out_dir + '/RightMedial_' + surf_name + '.png',
                       show=0)
