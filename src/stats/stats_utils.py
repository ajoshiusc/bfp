""" This module contains helpful utility function for running statistics using BFP """

import csv
import glob
import itertools
import multiprocessing
import os
import sys
from functools import partial

import numpy as np
import scipy as sp
import scipy.io as spio
import statsmodels.api as sm
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import fdrcorrection
from tqdm import tqdm

from brainsync import brainSync, normalizeData
from dfsio import readdfs, writedfs
from surfproc import (patch_color_attrib, smooth_patch, smooth_surf_function,
                      view_patch_vtk)

sys.path.append('../BrainSync')


def read_gord_data(data_dir, num_sub=1e6):

    dirlist = glob.glob(data_dir + '/*.filt.mat')
    subno = 0

    sub_data_files = []

    for fname in dirlist:
        full_fname = os.path.join(data_dir, fname)

        if os.path.isfile(full_fname) and subno < num_sub:
            sub_data_files.append(full_fname)
            subno += 1

    return sub_data_files


def read_demoCSV(csvfname, data_dir, file_ext, colsubj, colvar_exclude,
                 colvar_atlas, colvar_main, colvar_reg1, colvar_reg2):
    ''' loads csv file containing subjects' demographic information
        csv file should contain the following 5 columns for: subjectID, subjects to exclude (1=exclude), main effect variable, and 2 covariates to control for.
        if less than 2 covariates, create columns where all subjects have value of 1 so regression has no effect. '''
    file = open(csvfname)
    numline = len(file.readlines())
    subN = numline - 1

    sub_ID = []
    sub_fname = []
    subAtlas_idx = []
    reg_var = []
    reg_cvar1 = []
    reg_cvar2 = []
    count1 = 0
    pbar = tqdm(total=subN)
    lst = os.listdir(data_dir)
    with open(csvfname, newline='') as csvfile:
        creader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
        for row in creader:
            sub = row[colsubj]
            fname = os.path.join(data_dir, sub + "/func/" + sub + file_ext)
            if not os.path.isfile(fname) or int(row[colvar_exclude]) != 0:
                continue
            rvar = row[colvar_main]
            rcvar1 = row[colvar_reg1]
            rcvar2 = row[colvar_reg2]

            subAtlas_idx.append(row[colvar_atlas])
            sub_fname.append(fname)
            sub_ID.append(sub)
            reg_var.append(float(rvar))
            reg_cvar1.append(float(rcvar1))
            reg_cvar2.append(float(rcvar2))
            count1 += 1
            pbar.update(1)
            if count1 == subN:
                break

    pbar.close()
    print('CSV file read\nThere are %d subjects' % (len(sub_ID)))

    return sub_ID, sub_fname, subAtlas_idx, reg_var, reg_cvar1, reg_cvar2


def load_bfp_data(sub_fname, LenTime):
    ''' sub_fname: list of filenames of .mat files that contains Time x Vertex matrix of subjects' preprocessed fMRI data '''
    ''' LenTime: number of timepoints in data. this should be the same in all subjects '''
    ''' Outputs 3D matrix: Time x Vector x Subjects '''
    count1 = 0
    subN = len(sub_fname)
    print('loading data for ' + str(subN) + ' subjects')
    pbar = tqdm(total=subN)
    for ind in range(subN):
        fname = sub_fname[ind]
        df = spio.loadmat(fname)
        data = df['dtseries'].T
        if int(data.shape[0]) != LenTime:
            print(sub_fname[ind] +
                  ' does not have the correct number of timepoints')
        d, _, _ = normalizeData(data)

        if count1 == 0:
            sub_data = sp.zeros((LenTime, d.shape[1], subN))

        sub_data[:, :, count1] = d[:LenTime, ]
        count1 += 1
        pbar.update(1)
        if count1 == subN:
            break

    pbar.close()

    print('loaded data for ' + str(subN) + ' subjects')
    return sub_data

def sync2atlas(atlas, sub_data):
    print('Syncing to atlas, assume that the data is normalized')

    # Assume that the sub_data is already normalized
    syn_data = sp.zeros(sub_data.shape)
    for ind in tqdm(range(sub_data.shape[2])):
        syn_data[:, :, ind], _ = brainSync(X=atlas, Y=sub_data[:, :, ind])

    return syn_data

def dist2atlas(atlas, syn_data):
    ''' calculates geodesic distance between atlas and individual subjects at each vertex. all data should be synchronized to the atlas 
    inputs: atlas: Time x Vector matrix of reference atlas (see brainsync.py)
            syn_data: Time x Vector x Subjects matrix of subjects already synchronized to the atlas.
    output: diff Vector x Subjects data matrix'''
    numSub = syn_data.shape[2]
    numVert = syn_data.shape[1]
    print('calculating geodesic distances between ' + str(numSub) +
          ' subjects to the atlas in ' + str(numVert) + ' vertices.')
    count1 = 0
    pbar = tqdm(total=numSub)
    diff = sp.zeros([numVert, numSub])
    for ind in range(numSub):
        diff[:, ind] = sp.sum((syn_data[:, :, ind] - atlas)**2, axis=0)
        count1 += 1
        pbar.update(1)  # update the progress bar
        #print('%d,' % count1, end='')
        if count1 == numSub:
            break
    pbar.close()

    print('done')
    return diff


def sub2ctrl_dist(sub_file, ctrl_files, len_time=235):
    """ Compare a subject to controls """

    sub_data = spio.loadmat(sub_file)['dtseries'].T
    sub_data, _, _ = normalizeData(sub_data[:len_time, :])

    num_vert = sub_data.shape[1]
    fmri_diff = sp.zeros((num_vert, len(ctrl_files)))

    for ind, fname in enumerate(tqdm(ctrl_files)):
        ctrl_data = spio.loadmat(fname)['dtseries'].T
        ctrl_data, _, _ = normalizeData(ctrl_data[:len_time, :])
        ctrl_data, _ = brainSync(X=sub_data, Y=ctrl_data)
        fmri_diff[:, ind] = sp.sum((sub_data - ctrl_data)**2, axis=0)

    return fmri_diff


def pair_dist(rand_pair, sub_files, reg_var=[], len_time=235):
    """ Pair distance """
    sub1_data = spio.loadmat(sub_files[rand_pair[0]])['dtseries'].T
    sub2_data = spio.loadmat(sub_files[rand_pair[1]])['dtseries'].T

    sub1_data, _, _ = normalizeData(sub1_data[:len_time, :])
    sub2_data, _, _ = normalizeData(sub2_data[:len_time, :])

    sub2_data, _ = brainSync(X=sub1_data, Y=sub2_data)
    fmri_diff = sp.sum((sub2_data - sub1_data)**2, axis=0)
    if len(reg_var) > 0:
        regvar_diff = sp.square(reg_var[rand_pair[0]] - reg_var[rand_pair[1]])
        return fmri_diff, regvar_diff
    else:
        return fmri_diff


def pairsdist_regression(bfp_path,
                         sub_files,
                         reg_var,
                         num_perm=1000,
                         num_pairs=0,
                         len_time=235):
    """ Perform regression stats based on square distance between random pairs """

    # Get the number of vertices from a file
    num_vert = spio.loadmat(sub_files[0])['dtseries'].shape[0]
    num_sub = len(sub_files)

    # Allocate memory for subject data
    sub_data = np.zeros(shape=(len_time, num_vert, num_sub))

    #Generate random pairs
    print('Reading subjects')
    for subno, filename in enumerate(tqdm(sub_files)):
        data = spio.loadmat(filename)['dtseries'].T
        sub_data[:, :, subno], _, _ = normalizeData(data[:len_time, :])

    pairs = list(itertools.combinations(range(num_sub), r=2))

    if num_pairs > 0:
        rn = np.random.permutation(len(pairs))
        pairs = [pairs[i] for i in rn]
        pairs = pairs[:num_pairs]

    fmri_diff = sp.zeros((num_vert, len(pairs)))
    regvar_diff = sp.zeros(len(pairs))

    print('Computing pairwise differences')
    for pn, pair in enumerate(tqdm(pairs)):
        Y2, _ = brainSync(X=sub_data[:, :, pair[0]], Y=sub_data[:, :, pair[1]])
        fmri_diff[:, pn] = np.sum((Y2 - sub_data[:, :, pair[0]])**2, axis=0)
        regvar_diff[pn] = (reg_var[pair[0]] - reg_var[pair[1]])**2

    corr_pval = corr_perm_test(X=fmri_diff.T, Y=regvar_diff)

    #    corr_pval = sp.zeros(num_vert)
    #    for ind in tqdm(range(num_vert)):
    #        _, corr_pval[ind] = sp.stats.pearsonr(fmri_diff[ind, :], regvar_diff)
    #    corr_pval[sp.isnan(corr_pval)] = .5
    #

    labs = spio.loadmat(
        bfp_path +
        '/supp_data/USCBrain_grayordinate_labels.mat')['labels'].squeeze()
    labs[sp.isnan(labs)] = 0

    corr_pval[labs == 0] = 0.5

    corr_pval_fdr = 0.5 * sp.ones(num_vert)
    _, corr_pval_fdr[labs > 0] = fdrcorrection(corr_pval[labs > 0])

    return corr_pval, corr_pval_fdr


def corr_pearson_fdr(X, Y, nperm=1000):
    #X: nsub x vertices
    #Y: cognitive scores nsub X 1
    num_vert = X.shape[1]

    corr_pval = sp.zeros(num_vert)
    for ind in tqdm(range(num_vert)):
        _, corr_pval[ind] = sp.stats.pearsonr(X[:, ind], Y)

    corr_pval[sp.isnan(corr_pval)] = .5

    _, corr_pval_fdr = fdrcorrection(corr_pval)

    return corr_pval_fdr


def corr_perm_test(X, Y, nperm=1000):
    #X: nsub x vertices
    #Y: cognitive scores nsub X 1

    X, _, _ = normalizeData(X)
    Y, _, _ = normalizeData(Y)

    nsub = X.shape[0]
    rho_vert = np.sum(X * Y[:, None], axis=0)
    max_null = np.zeros(nperm)

    print('Permutation testing')
    for ind in tqdm(range(nperm)):
        perm1 = np.random.permutation(nsub)
        max_null[ind] = np.amax(np.sum(X * Y[perm1, None], axis=0))

    pval = np.sum(rho_vert[:, None] < max_null[None, :], axis=1) / nperm

    return pval


def randpairsdist_reg_parallel(bfp_path,
                               sub_files,
                               reg_var,
                               num_pairs=2000,
                               nperm=1000,
                               len_time=235,
                               num_proc=4,
                               fdr_test=False):
    """ Perform regression stats based on square distance between random pairs """

    # Get the number of vertices from a file
    num_vert = spio.loadmat(sub_files[0])['dtseries'].shape[0]

    # Generate pairs
    pairs = list(itertools.combinations(range(len(sub_files)), r=2))

    if num_pairs > 0:
        rn = np.random.permutation(len(pairs))
        pairs = [pairs[i] for i in rn]
        if num_pairs < len(pairs):
            pairs = pairs[:num_pairs]
        else:
            num_pairs = len(pairs)

    fmri_diff = sp.zeros((num_vert, num_pairs))
    regvar_diff = sp.zeros(num_pairs)

    results = multiprocessing.Pool(num_proc).imap(
        partial(
            pair_dist, sub_files=sub_files, reg_var=reg_var,
            len_time=len_time), pairs)

    ind = 0
    for res in results:
        fmri_diff[:, ind] = res[0]
        regvar_diff[ind] = res[1]
        ind += 1

    if not fdr_test:
        print('Performing Permutation test with MAX statistic')
        corr_pval = corr_perm_test(X=fmri_diff.T, Y=regvar_diff, nperm=nperm)
    else:
        print('Performing Pearson correlation with FDR testing')
        corr_pval = corr_pearson_fdr(X=fmri_diff.T, Y=regvar_diff, nperm=nperm)

    corr_pval[sp.isnan(corr_pval)] = .5

    labs = spio.loadmat(
        bfp_path +
        '/supp_data/USCBrain_grayordinate_labels.mat')['labels'].squeeze()
    labs[sp.isnan(labs)] = 0

    corr_pval[labs == 0] = 0.5

    return corr_pval


def group_diff_fdr(grp1, grp2, alt_hypo='less'):
    '''Input grp1 = numvert x (numsamples grp1), grp2 = numvert x (numsamples grp2)
    alt_hypo: 'less', 'more' or 'two-sided'''

    print('Performing Mann Whitney test, checking grp1<grp2')
    pval = sp.zeros(grp1.shape[0])

    for vind in tqdm(range(grp1.shape[0])):

        _, pval[vind] = sp.stats.ranksums(grp1[vind, :], grp2[vind, :])


#        _, pval[vind] = sp.stats.mannwhitneyu(
#            grp1[vind, :],
#            grp2[vind, :],
#            use_continuity=True,
#            alternative=alt_hypo)

    _, pval_fdr = fdrcorrection(pval)
    return pval_fdr, pval


def compare_sub2ctrl(bfp_path,
                     sub_file,
                     ctrl_files,
                     num_pairs=2000,
                     nperm=1000,
                     len_time=235,
                     num_proc=4,
                     fdr_test=True):

    # Get the number of vertices from a file
    num_vert = spio.loadmat(ctrl_files[0])['dtseries'].shape[0]

    # Generate pairs
    pairs = list(itertools.combinations(range(len(ctrl_files)), r=2))

    if num_pairs > 0:
        rn = np.random.permutation(len(pairs))
        pairs = [pairs[i] for i in rn]
        if num_pairs < len(pairs):
            pairs = pairs[:num_pairs]
        else:
            num_pairs = len(pairs)

    fmri_diff_null = sp.zeros((num_vert, num_pairs))

    if num_proc == 1:
        for ind in tqdm(range(len(pairs))):
            fmri_diff_null[:, ind] = pair_dist(
                sub_files=ctrl_files, len_time=len_time, rand_pair=pairs[ind])

    else:
        results = multiprocessing.Pool(num_proc).imap(
            partial(pair_dist, sub_files=ctrl_files, len_time=len_time), pairs)

        ind = 0
        for res in results:
            fmri_diff_null[:, ind] = res[0]
            ind += 1

    sub2ctrl_diff = sub2ctrl_dist(
        sub_file=sub_file, ctrl_files=ctrl_files, len_time=len_time)

    if not fdr_test:
        print('Performing Permutation test with MAX statistic')
        #corr_pval = corr_perm_test(X=fmri_diff.T, Y=[], nperm=nperm)
    else:
        print('Performing Pearson correlation with FDR testing')
        pval_fdr, pval = group_diff_fdr(
            grp1=fmri_diff_null, grp2=sub2ctrl_diff, alt_hypo='less')

    return pval_fdr, pval


'''Deprecated'''


def randpairsdist_reg(bfp_path,
                      sub_files,
                      reg_var,
                      num_pairs=1000,
                      len_time=235):
    """ Perform regression stats based on square distance between random pairs """
    print('dist2atlas_reg, assume that the data is normalized')
    print('This function is deprecated!!!!!!!!!!')

    # Get the number of vertices from a file
    num_vert = spio.loadmat(sub_files[0])['dtseries'].shape[0]

    #Generate random pairs
    rand_pairs = sp.random.choice(len(sub_files), (num_pairs, 2), replace=True)

    fmri_diff = sp.zeros((num_vert, num_pairs))
    regvar_diff = sp.zeros(num_pairs)

    print('Reading subjects')

    # Compute distance to atlas
    for ind in tqdm(range(num_pairs)):
        sub1_data = spio.loadmat(sub_files[rand_pairs[ind, 0]])['dtseries'].T
        sub2_data = spio.loadmat(sub_files[rand_pairs[ind, 1]])['dtseries'].T

        sub1_data, _, _ = normalizeData(sub1_data[:len_time, :])
        sub2_data, _, _ = normalizeData(sub2_data[:len_time, :])

        sub2_data, _ = brainSync(X=sub1_data, Y=sub2_data)
        fmri_diff[:, ind] = sp.sum((sub2_data - sub1_data)**2, axis=0)
        regvar_diff[ind] = sp.square(reg_var[rand_pairs[ind, 0]] -
                                     reg_var[rand_pairs[ind, 1]])

    corr_pval = sp.zeros(num_vert)
    for ind in tqdm(range(num_vert)):
        _, corr_pval[ind] = sp.stats.pearsonr(fmri_diff[ind, :], regvar_diff)

    corr_pval[sp.isnan(corr_pval)] = .5

    labs = spio.loadmat(bfp_path + '/supp_data/USCBrain_grayord_labels.mat'
                        )['labels'].squeeze()

    corr_pval_fdr = sp.zeros(num_vert)
    _, corr_pval_fdr[labs > 0] = fdrcorrection(corr_pval[labs > 0])

    return corr_pval, corr_pval_fdr


def sync2atlas(atlas, sub_data):
    print('Syncing to atlas, assume that the data is normalized')

    # Assume that the sub_data is already normalized
    syn_data = sp.zeros(sub_data.shape)
    for ind in tqdm(range(sub_data.shape[2])):
        syn_data[:, :, ind], _ = brainSync(X=atlas, Y=sub_data[:, :, ind])

    return syn_data


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


def LinReg_resid(x, y):
    slope, intercept, _, _, _ = sp.stats.linregress(x, y)
    predicted = x * slope + intercept
    resid = y - predicted

    return resid

def multiLinReg_resid(x,y):
    regr = sklearn.linear_model.LinearRegression()
    regr.fit(x,y)
    resid = y - regr.predict(x)
    return resid

def multiLinReg_corr(subTest_diff, subTest_varmain, subTest_varc1, subTest_varc2):
    subTest_varc12 = sp.zeros((subTest_varc1.shape[0],2))
    for i in range(subTest_varc1.shape[0]):
        subTest_varc12[i,0] = subTest_varc1[i]
        subTest_varc12[i,1] = subTest_varc2[i]
    print('regressing out 2 covariates')
    diff_resid1 = sp.zeros(subTest_diff.shape)
    numV = subTest_diff.shape[0]
    for nv in tqdm(range(numV)):
        diff_resid1[nv, :] = multiLinReg_resid(subTest_varc12, subTest_diff[nv, :])
        
    print('computing correlation against main variable')
    rval = sp.zeros(numV)
    pval = sp.zeros(numV)
    for nv in tqdm(range(numV)):
        rval[nv], pval[nv] = sp.stats.pearsonr(
            subTest_varmain, diff_resid1[nv, :])

    a = spio.loadmat('supp_data/USCBrain_grayordinate_labels.mat')
    labs = a['labels'].squeeze()
    labs[sp.isnan(labs)] = 0
    pval_fdr = sp.zeros(numV)
    _, pv = fdrcorrection(pval[labs > 0])
    pval_fdr[labs > 0] = pv

    return rval, pval, pval_fdr

def LinReg_corr(subTest_diff, subTest_varmain, subTest_varc1, subTest_varc2):
    print('regressing out 1st covariate')
    diff_resid1 = sp.zeros(subTest_diff.shape)
    numV = subTest_diff.shape[0]
    for nv in tqdm(range(numV)):
        diff_resid1[nv, :] = LinReg_resid(subTest_varc1, subTest_diff[nv, :])

    print('regressing out 2nd covariate')
    diff_resid2 = sp.zeros(subTest_diff.shape)
    for nv in tqdm(range(numV)):
        diff_resid2[nv, :] = LinReg_resid(subTest_varc2, diff_resid1[nv, :])

    print('computing correlation against main variable')
    rval = sp.zeros(numV)
    pval = sp.zeros(numV)
    for nv in tqdm(range(numV)):
        _, _, rval[nv], pval[nv], _ = sp.stats.linregress(
            subTest_varmain, diff_resid2[nv, :])

    a = spio.loadmat('supp_data/USCBrain_grayordinate_labels.mat')
    labs = a['labels'].squeeze()
    labs[sp.isnan(labs)] = 0
    pval_fdr = sp.zeros(numV)
    _, pv = fdrcorrection(pval[labs > 0])
    pval_fdr[labs > 0] = pv

    return rval, pval, pval_fdr