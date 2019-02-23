""" This module contains helpful utility function for running statistics using BFP """

import csv
import os
import scipy as sp
import scipy.io as spio
from tqdm import tqdm
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
from sklearn.decomposition import PCA
from surfproc import view_patch_vtk, patch_color_attrib, smooth_surf_function, smooth_patch
from dfsio import readdfs
from brainsync import normalizeData, brainSync


def read_fcon1000_data(csv_fname,
                       data_dir,
                       reg_var_name='Verbal IQ',
                       num_sub=5,
                       reg_var_positive=1):
    """ reads fcon1000 csv and data"""

    count1 = 0
    sub_ids = []
    reg_var = []
    pbar = tqdm(total=num_sub)

    with open(csv_fname, newline='') as csvfile:
        creader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
        for row in creader:

            # read the regression variable
            rvar = row[reg_var_name]

            # Read the filtered data by default
            fname = os.path.join(
                data_dir, row['ScanDir ID'] + '_rest_bold.32k.GOrd.filt.mat')

            # If the data does not exist for this subject then skip it
            if not os.path.isfile(fname) or int(row['QC_Rest_1']) != 1:
                continue

            if reg_var_positive == 1 and sp.float64(rvar) < 0:
                continue

            if count1 == 0:
                sub_data_files = []

            # Truncate the data at a given number of time samples This is needed because
            # BrainSync needs same number of time sampples
            sub_data_files.append(fname)
            sub_ids.append(row['ScanDir ID'])
            reg_var.append(float(rvar))

            count1 += 1
            pbar.update(1)  # update the progress bar
            #print('%d,' % count1, end='')
            if count1 == num_sub:
                break

    pbar.close()
    print('CSV file and the data has been read\nThere are %d subjects' %
          (len(sub_ids)))

    return sub_ids, sp.array(reg_var), sub_data_files


def sync2atlas(atlas, sub_data):
    print('Syncing to atlas, assume that the data is normalized')

    # Assume that the sub_data is already normalized
    syn_data = sp.zeros(sub_data.shape)
    for ind in tqdm(range(sub_data.shape[2])):
        syn_data[:, :, ind], _ = brainSync(X=atlas, Y=sub_data[:, :, ind])

    return syn_data


def ref_avg_atlas(ref_id, sub_files, len_time=235):
    ''' Generates atlas by syncing to one reference subject'''

    ref_data = spio.loadmat(sub_files[ref_id])['dtseries'].T
    ref_data, _, _ = normalizeData(ref_data[:len_time, :])

    for ind in tqdm(range(len(sub_files))):
        sub_data = spio.loadmat(sub_files[ind])['dtseries'].T
        sub_data, _, _ = normalizeData(sub_data[:len_time, :])
        s_data, _ = brainSync(X=ref_data, Y=sub_data)
        if ind == 0:
            avg_atlas = s_data
        else:
            avg_atlas += s_data

    avg_atlas /= len(sub_files)

    #    avg_atlas, _, _ = normalizeData(avg_atlas)

    return avg_atlas


def dist2atlas_reg(bfp_path, ref_atlas, sub_files, reg_var, len_time=235):
    """ Perform regression stats based on square distance to atlas """
    print('dist2atlas_reg, assume that the data is normalized')

    num_vert = ref_atlas.shape[1]
    num_sub = len(sub_files)

    # Take absolute value of difference from the mean
    # for the IQ measure
    reg_var = sp.absolute(reg_var - sp.mean(reg_var))

    diff = sp.zeros((num_vert, num_sub))

    # Compute distance to atlas
    for ind in tqdm(range(num_sub)):
        sub_data = spio.loadmat(sub_files[ind])['dtseries'].T
        sub_data, _, _ = normalizeData(sub_data[:len_time, :])
        Y2, _ = brainSync(X=ref_atlas, Y=sub_data)
        diff[:, ind] = sp.sum((Y2 - ref_atlas)**2, axis=0)

    corr_pval = sp.zeros(num_vert)
    for vrt in tqdm(range(num_vert)):
        _, corr_pval[vrt] = sp.stats.pearsonr(diff[vrt, :], reg_var)

    corr_pval[sp.isnan(corr_pval)] = .5

    lab = spio.loadmat(bfp_path + '/supp_data/USCBrain_grayord_labels.mat')
    labs = lab['labels'].squeeze()

    corr_pval_fdr = sp.zeros(num_vert)
    _, pv = fdrcorrection(corr_pval[labs > 0])
    corr_pval_fdr[labs > 0] = pv

    return corr_pval, corr_pval_fdr


def lin_reg(bfp_path,
            ref_atlas,
            sub_files,
            reg_var,
            Vndim=235,
            Sndim=20,
            len_time=235):
    """ Perform regression stats based on distance to atlas """

    num_vert = ref_atlas.shape[1]
    num_sub = len(sub_files)
    a = spio.loadmat(bfp_path + '/supp_data/USCBrain_grayord_labels.mat')
    labs = a['labels'].squeeze()

    labs[sp.isnan(labs)] = 0
    print('Computing PCA basis function from the atlas')
    pca = PCA(n_components=Vndim)
    pca.fit(ref_atlas.T)

    reduced_data = sp.zeros((Vndim, num_vert, num_sub))
    for ind in tqdm(range(num_sub)):

        sub_data = spio.loadmat(sub_files[ind])['dtseries'].T
        sub_data, _, _ = normalizeData(sub_data[:len_time, :])
        Y2, _ = brainSync(X=ref_atlas, Y=sub_data)

        if Vndim == len_time:
            reduced_data[:, :, ind] = sub_data
        else:
            reduced_data[:, :, ind] = pca.transform(Y2.T).T

    pval_linreg = sp.zeros(num_vert)

    pca = PCA(n_components=Sndim)

    for vrt in tqdm(range(num_vert)):
        X = reduced_data[:, vrt, :]
        if Sndim != num_sub:
            pca.fit(X.T)
            X = pca.transform(X.T).T
        X = sm.add_constant(X.T)
        est = sm.OLS(reg_var, X)
        pval_linreg[vrt] = est.fit().f_pvalue

    print('Regression is done')

    pval_linreg[sp.isnan(pval_linreg)] = .5

    pval_linreg_fdr = sp.zeros(num_vert)
    _, pv = fdrcorrection(pval_linreg[labs > 0])
    pval_linreg_fdr[labs > 0] = pv

    return pval_linreg, pval_linreg_fdr


def vis_save_pval(bfp_path, pval_map, surf_name, smooth_iter=1500):
    lsurf = readdfs(bfp_path + '/supp_data/bci32kleft.dfs')
    rsurf = readdfs(bfp_path + '/supp_data/bci32kright.dfs')

    lsurf.attributes = sp.zeros((lsurf.vertices.shape[0]))
    rsurf.attributes = sp.zeros((rsurf.vertices.shape[0]))
    lsurf = smooth_patch(lsurf, iterations=smooth_iter)
    rsurf = smooth_patch(rsurf, iterations=smooth_iter)

    num_vert = lsurf.vertices.shape[0]

    lsurf.attributes = 0.05 - pval_map.squeeze()
    lsurf.attributes = lsurf.attributes[:num_vert]
    rsurf.attributes = 0.05 - pval_map.squeeze()
    rsurf.attributes = rsurf.attributes[num_vert:2 * num_vert]

    lsurf = patch_color_attrib(lsurf, clim=[0, .05])
    rsurf = patch_color_attrib(rsurf, clim=[0, .05])

    # If p value above .05 then make the surface grey
    lsurf.vColor[lsurf.attributes < 0, :] = .5
    rsurf.vColor[rsurf.attributes < 0, :] = .5

    # Visualize left hemisphere
    view_patch_vtk(
        lsurf,
        azimuth=100,
        elevation=180,
        roll=90,
        outfile='left_' + surf_name + '_pval.png',
        show=0)
    # Visualize right hemisphere
    view_patch_vtk(
        rsurf,
        azimuth=-100,
        elevation=180,
        roll=-90,
        outfile='right_' + surf_name + '_pval.png',
        show=0)
