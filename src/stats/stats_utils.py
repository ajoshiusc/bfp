""" This module contains helpful utility function for running statistics using BFP """

# try:
#    import vtk
#    VTK_INSTALLED = 1
# except ImportError as e:
#    VTK_INSTALLED = 0
#VTK_INSTALLED = 0

import csv
import glob
import itertools
import multiprocessing
import os
import sys
from functools import partial
from itertools import product
from multiprocessing import Pool

import numpy as np
import scipy as sp
import scipy.io as spio
import scipy.stats as ss
import sklearn
import statsmodels.api as sm
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import fdrcorrection
from tqdm import tqdm

from brainsync import brainSync, normalizeData
from dfsio import readdfs, writedfs
from surfproc import patch_color_attrib, smooth_surf_function

# if VTK_INSTALLED:
#    from surfproc import view_patch_vtk, smooth_patch


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

def read_bord_data(data_dir, num_sub=1e6):

    dirlist = glob.glob(data_dir + '/*BOrd.mat')
    subno = 0

    sub_data_files = []

    for fname in dirlist:
        full_fname = os.path.join(data_dir, fname)

        if os.path.isfile(full_fname) and subno < num_sub:
            sub_data_files.append(full_fname)
            subno += 1

    return sub_data_files

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
    output: geo_dist Vector x Subjects data matrix of Geodesic distances
            pearson_corr Vector x Subjects data matrix of Pearson correlations
    '''

    numSub = syn_data.shape[2]
    numVert = syn_data.shape[1]
    print('calculating geodesic distances between ' + str(numSub) +
          ' subjects to the atlas in ' + str(numVert) + ' vertices.')
    count1 = 0
    pbar = tqdm(total=numSub)
    pearson_corr = sp.zeros([numVert, numSub])
    geo_dist = sp.zeros([numVert, numSub])

    for ind in range(numSub):
        pearson_corr[:, ind] = sp.sum((syn_data[:, :, ind] * atlas), axis=0)
        geo_dist[:, ind] = np.arccos(pearson_corr[:, ind])

        count1 += 1
        pbar.update(1)  # update the progress bar
        #print('%d,' % count1, end='')
        if count1 == numSub:
            break
    pbar.close()
    a = pearson_corr > 1
    geo_dist[a] = 0
    print('done')
    return geo_dist, pearson_corr


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


def pair_dist_two_groups(rand_pair,
                         sub_grp1_files,
                         sub_grp2_files,
                         sub_data1=[],
                         sub_data2=[],
                         len_time=235):

    sub_data1 = np.array(sub_data1)
    sub_data2 = np.array(sub_data2)
    """ Pair distance for two groups of subjects """
    if sub_data1.size > 0:
        sub1_data = sub_data1[:, :, rand_pair[0]]
    else:
        sub1_data = spio.loadmat(sub_grp1_files[rand_pair[0]])['dtseries'].T
        sub1_data, _, _ = normalizeData(sub1_data[:len_time, :])

    if sub_data2.size > 0:
        sub2_data = sub_data2[:, :, rand_pair[1]]
    else:
        sub2_data = spio.loadmat(sub_grp2_files[rand_pair[1]])['dtseries'].T
        sub2_data, _, _ = normalizeData(sub1_data[:len_time, :])

    sub2_data, _ = brainSync(X=sub1_data, Y=sub2_data)
    fmri_diff = sp.sum((sub2_data - sub1_data)**2, axis=0)

    # Returns SQUARE of the distance
    return fmri_diff


def pair_dist(rand_pair, sub_files, sub_data=[], reg_var=[], len_time=235):
    """ Pair distance """
    sub_data = np.array(sub_data)
    if sub_data.size > 0:
        sub1_data = sub_data[:, :, rand_pair[0]]
        sub2_data = sub_data[:, :, rand_pair[1]]
    else:
        sub1_data = spio.loadmat(sub_files[rand_pair[0]])['dtseries'].T
        sub2_data = spio.loadmat(sub_files[rand_pair[1]])['dtseries'].T
        sub1_data, _, _ = normalizeData(sub1_data[:len_time, :])
        sub2_data, _, _ = normalizeData(sub2_data[:len_time, :])

    sub2_data, _ = brainSync(X=sub1_data, Y=sub2_data)
    fmri_diff = sp.sum((sub2_data - sub1_data)**2, axis=0)

    # Returns SQUARE of the distance
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

    # Generate random pairs
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


def corr_pearson_fdr(X_pairs, Y_pairs, reg_var, num_sub, nperm=1000):
    # X: nsub x vertices
    # Y: cognitive scores nsub X 1
    X, _, _ = normalizeData(X_pairs)

    Y, _, _ = normalizeData(Y_pairs[:, None])

    num_vert = X.shape[1]

    corr_pval = sp.zeros(num_vert)
    for ind in tqdm(range(num_vert)):
        _, corr_pval[ind] = sp.stats.pearsonr(X[:, ind], Y.squeeze())

    corr_pval[sp.isnan(corr_pval)] = .5

    _, corr_pval_fdr = fdrcorrection(corr_pval)

    return corr_pval_fdr, corr_pval


def corr_perm_test(X_pairs, Y_pairs, reg_var, num_sub, nperm=1000):
    # X: nsub x vertices
    # Y: cognitive scores nsub X 1

    X, _, _ = normalizeData(X_pairs)

    num_pairs = X.shape[0]
    Y_pairs, _, _ = normalizeData(Y_pairs[:, None])
    rho_orig = np.sum(X * Y_pairs, axis=0)
    max_null = np.zeros(nperm)
    n_count = np.zeros(X.shape[1])

    print('Permutation testing')
    for ind in tqdm(range(nperm)):
        pairs, _ = gen_rand_pairs(num_sub=num_sub, num_pairs=num_pairs)
        pairs = np.array(pairs)
        Y = sp.square(reg_var[pairs[:, 0]] - reg_var[pairs[:, 1]])

        Y, _, _ = normalizeData(Y[:, None])

        rho_perm = np.sum(X * Y, axis=0)
        max_null[ind] = np.amax(rho_perm)
        n_count += np.float32(rho_perm >= rho_orig)

    pval_max = np.sum(rho_orig[:, None] <= max_null[None, :], axis=1) / nperm

    pval_perm = n_count / nperm

    _, pval_perm_fdr = fdrcorrection(pval_perm)

    return pval_max, pval_perm_fdr, pval_perm


'''
#
def rand_pair_dist(sub_files=sub_fname, num_pairs=num_pairs):

    print('hi')

    num_sub = len(sub_fname)
    pairs = gen_rand_pairs(num_sub, num_pairs)

    return pairs, pair_dist

'''


def gen_rand_pairs(num_sub, num_pairs):
    # Generate pairs
    pairs = list(itertools.combinations(range(num_sub), r=2))

    if num_pairs > 0:
        rn = np.random.permutation(len(pairs))
        pairs = [pairs[i] for i in rn]
        if num_pairs < len(pairs):
            pairs = pairs[:num_pairs]
        else:
            num_pairs = len(pairs)

    return pairs, num_pairs


def randpair_groupdiff(sub_grp1_files, sub_grp2_files, num_pairs,
                       len_time=255):

    print('Grp diff')

    num_vert = spio.loadmat(sub_grp1_files[0])['dtseries'].shape[0]

    print('Generating random pairs from group 1')
    pairs_grp1, num_pairs1 = gen_rand_pairs(num_sub=len(sub_grp1_files),
                                            num_pairs=num_pairs)

    fmri_diff1 = sp.zeros((num_vert, num_pairs1))

    # Preload data This only slighly faster, better is to load on the fly and multiprocess
    print('Reading data for group 1')
    sub_data1 = np.zeros((len_time, num_vert, len(sub_grp1_files)))
    for i, fname in enumerate(tqdm(sub_grp1_files)):
        sub1_data = spio.loadmat(fname)['dtseries'][:, :len_time].T
        sub_data1[:, :, i], _, _ = normalizeData(sub1_data)

    print('Compute differences in fMRI of random pairs from group 1')
    for i, rand_pair in enumerate(tqdm(pairs_grp1)):
        fmri_diff1[:, i] = pair_dist(rand_pair=rand_pair,
                                     sub_files=sub_grp1_files,
                                     sub_data=sub_data1,
                                     len_time=len_time)

    S1 = 0.5 * np.mean(fmri_diff1, axis=1)

    print('Generating random pairs from group 2')
    pairs_grp2, num_pairs2 = gen_rand_pairs(num_sub=len(sub_grp2_files),
                                            num_pairs=num_pairs)

    fmri_diff2 = sp.zeros((num_vert, num_pairs2))

    # Preload data for group 2
    print('Reading data for group 2')
    sub_data2 = np.zeros((len_time, num_vert, len(sub_grp2_files)))
    for i, fname in enumerate(tqdm(sub_grp2_files)):
        sub2_data = spio.loadmat(fname)['dtseries'][:, :len_time].T
        sub_data2[:, :, i], _, _ = normalizeData(sub2_data)

    print('Compute differences in fMRI of random pairs from group 2')
    for i, rand_pair in enumerate(tqdm(pairs_grp2)):
        fmri_diff2[:, i] = pair_dist(rand_pair=rand_pair,
                                     sub_files=sub_grp2_files,
                                     sub_data=sub_data2,
                                     len_time=len_time)

    S2 = 0.5 * np.mean(fmri_diff2, axis=1)

    print('Generating random pairs from all subjects (grp1 + grp2)')

    # Generating random pairs. For large group this may allocate huge amount of memory,
    # use the following solution in that case
    # https://stackoverflow.com/questions/36779729/shuffling-combinations-without-converting-iterable-itertools-combinations-to-l

    all_pairs = np.array(
        list(product(range(len(sub_grp1_files)), range(len(sub_grp2_files)))))
    sp.random.shuffle(all_pairs)
    all_pairs = all_pairs[:num_pairs, :]
    fmri_diff = sp.zeros((num_vert, all_pairs.shape[0]))

    print('Compute differences in fMRI of random pairs from group1 to group 2')
    for i, rand_pair in enumerate(tqdm(all_pairs)):
        fmri_diff[:, i] = pair_dist_two_groups(rand_pair=rand_pair,
                                               sub_grp1_files=sub_grp1_files,
                                               sub_grp2_files=sub_grp2_files,
                                               sub_data1=sub_data1,
                                               sub_data2=sub_data2,
                                               len_time=len_time)

    # We will perform Welch's t test (modified in a pairwise stats)
    # https://en.wikipedia.org/wiki/Welch%27s_t-test

    n1 = sub_data1.shape[2]
    n2 = sub_data2.shape[2]

    tscore = np.sqrt((0.5 * np.mean(fmri_diff, axis=1) - (S1 * n1 + S2 * n2) /
                      (n1 + n2))) / np.sqrt(S1 / n1 + S1 / n2 + 1e-6)

    tscore[np.isnan(tscore)] = 0

    dof = (S1 / n1 + S2 / n2 + 1e-6)**2 / (S1**2 / ((n1**2) * (n1 - 1)) + S2**2 /
                                           ((n2**2) * (n2 - 1)) + 1e-6)
    pval = sp.stats.t.sf(tscore, dof) * 2  # two-sided pvalue = Prob(abs(t)>tt)

    return tscore, pval


def randpair_groupdiff_ftest(sub_grp1_files, sub_grp2_files, num_pairs,
                             len_time=255):

    print('Grp diff using f-test and brainsync')

    num_vert = spio.loadmat(sub_grp1_files[0])['dtseries'].shape[0]

    print('Generating random pairs from group 1')
    pairs_grp1, num_pairs1 = gen_rand_pairs(num_sub=len(sub_grp1_files),
                                            num_pairs=num_pairs)

    fmri_diff1 = sp.zeros((num_vert, num_pairs1))

    # Preload data This only slighly faster, better is to load on the fly and multiprocess
    print('Reading data for group 1')
    sub_data1 = np.zeros((len_time, num_vert, len(sub_grp1_files)))
    for i, fname in enumerate(tqdm(sub_grp1_files)):
        sub1_data = spio.loadmat(fname)['dtseries'][:, :len_time].T
        sub_data1[:, :, i], _, _ = normalizeData(sub1_data)

    print('Compute differences in fMRI of random pairs from group 1')
    for i, rand_pair in enumerate(tqdm(pairs_grp1)):
        fmri_diff1[:, i] = pair_dist(rand_pair=rand_pair,
                                     sub_files=sub_grp1_files,
                                     sub_data=sub_data1,
                                     len_time=len_time)

    S1 = 0.5 * np.mean(fmri_diff1, axis=1)

    print('Generating random pairs from group 2')
    pairs_grp2, num_pairs2 = gen_rand_pairs(num_sub=len(sub_grp2_files),
                                            num_pairs=num_pairs)

    fmri_diff2 = sp.zeros((num_vert, num_pairs2))

    # Preload data for group 2
    print('Reading data for group 2')
    sub_data2 = np.zeros((len_time, num_vert, len(sub_grp2_files)))
    for i, fname in enumerate(tqdm(sub_grp2_files)):
        sub2_data = spio.loadmat(fname)['dtseries'][:, :len_time].T
        sub_data2[:, :, i], _, _ = normalizeData(sub2_data)

    print('Compute differences in fMRI of random pairs from group 2')
    for i, rand_pair in enumerate(tqdm(pairs_grp2)):
        fmri_diff2[:, i] = pair_dist(rand_pair=rand_pair,
                                     sub_files=sub_grp2_files,
                                     sub_data=sub_data2,
                                     len_time=len_time)

    S2 = 0.5 * np.mean(fmri_diff2, axis=1)

    # We will perform f-test test (modified in a pairwise stats)
    #

    n1 = sub_data1.shape[2]*len_time
    n2 = sub_data2.shape[2]*len_time

    F = S1 / (S2 + 1e-16)

    pval = 1-ss.f.cdf(F, n1 - 1, n2 - 1)

    return F, pval


def randpairs_regression(bfp_path,
                         sub_files,
                         reg_var,
                         num_pairs=2000,
                         nperm=1000,
                         len_time=235,
                         num_proc=4,
                         pearson_fdr_test=False):
    """ Perform regression stats based on square distance between random pairs """

    # Get the number of vertices from a file
    num_vert = spio.loadmat(sub_files[0])['dtseries'].shape[0]

    pairs, num_pairs = gen_rand_pairs(num_sub=len(sub_files),
                                      num_pairs=num_pairs)

    fmri_diff = sp.zeros((num_vert, num_pairs))
    regvar_diff = sp.zeros(num_pairs)

    if num_proc > 1:
        pool = Pool(num_proc)

        results = pool.imap(
            partial(pair_dist,
                    sub_files=sub_files,
                    reg_var=reg_var,
                    len_time=len_time), pairs)

        ind = 0
        for res in results:
            fmri_diff[:, ind] = res[0]
            regvar_diff[ind] = res[1]
            ind += 1

    else:
        for ind in tqdm(range(len(pairs))):

            fmri_diff[:, ind], regvar_diff[ind] = pair_dist(
                sub_files=sub_files,
                reg_var=reg_var,
                len_time=len_time,
                rand_pair=pairs[ind])

    corr_pval2 = 0
    if not pearson_fdr_test:
        print('Performing Permutation test with MAX statistic')
        corr_pval, corr_pval2, _ = corr_perm_test(X_pairs=fmri_diff.T,
                                                  Y_pairs=regvar_diff,
                                                  reg_var=reg_var,
                                                  num_sub=len(sub_files),
                                                  nperm=nperm)
    else:
        print('Performing Pearson correlation with FDR testing')
        corr_pval, corr_pval2 = corr_pearson_fdr(X_pairs=fmri_diff.T,
                                                 Y_pairs=regvar_diff,
                                                 reg_var=reg_var,
                                                 num_sub=len(sub_files),
                                                 nperm=nperm)

    corr_pval[sp.isnan(corr_pval)] = .5

    labs = spio.loadmat(
        bfp_path +
        '/supp_data/USCBrain_grayordinate_labels.mat')['labels'].squeeze()
    labs[sp.isnan(labs)] = 0

    corr_pval[labs == 0] = 0.5

    return corr_pval, corr_pval2


def dist_diff_fdr(grp1, grp2):
    '''Input grp1 = numvert x (numsamples grp1), grp2 = numvert x (numsamples grp2)
    alt_hypo: 'less', 'more' or 'two-sided'''  # These are distances so need to use permutation test
    # Tests if Alt Hypo grp2 > grp1

    pval = 0
    for j in range(grp2.shape[1]):
        pval += (grp1[:, :] > grp2[:, j][:, None]).sum(axis=1)/grp1.shape[1]

    pval /= grp2.shape[1]

    _, pval_fdr = fdrcorrection(pval)

    return pval_fdr, pval


def group_diff_fdr(grp1, grp2, alt_hypo='less'):
    '''Input grp1 = numvert x (numsamples grp1), grp2 = numvert x (numsamples grp2)
    alt_hypo: 'less', 'more' or 'two-sided'''

    print('Performing Mann Whitney test, checking grp1<grp2')
    pval = sp.zeros(grp1.shape[0])

    for vind in tqdm(range(grp1.shape[0])):

        #        _, pval[vind] = sp.stats.ranksums(grp1[vind, :], grp2[vind, :])

        if np.linalg.norm(grp2[vind, :]) > 0:

            _, pval[vind] = sp.stats.mannwhitneyu(
                grp1[vind, :],
                grp2[vind, :],
                use_continuity=True,
                alternative=alt_hypo)
        else:
            pval[vind] = 0.5

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
            fmri_diff_null[:, ind] = pair_dist(sub_files=ctrl_files,
                                               len_time=len_time,
                                               rand_pair=pairs[ind])

    else:
        results = multiprocessing.Pool(num_proc).imap(
            partial(pair_dist, sub_files=ctrl_files, len_time=len_time), pairs)

        ind = 0
        for res in results:
            fmri_diff_null[:, ind] = res[0]
            ind += 1

    sub2ctrl_diff = sub2ctrl_dist(sub_file=sub_file,
                                  ctrl_files=ctrl_files,
                                  len_time=len_time)

    if not fdr_test:
        print('Performing Permutation test with MAX statistic')
        #corr_pval = corr_perm_test(X=fmri_diff.T, Y=[], nperm=nperm)
    else:
        print('Performing Pearson correlation with FDR testing')
        pval_fdr, pval = group_diff_fdr(grp1=fmri_diff_null,
                                        grp2=sub2ctrl_diff)
#                                        alt_hypo='less')

    return pval_fdr, pval


'''Deprecated'''


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


def multiLinReg_resid(x, y):
    regr = sklearn.linear_model.LinearRegression()
    regr.fit(x, y)
    resid = y - regr.predict(x)
    return resid


def multiLinReg_corr(subTest_diff, subTest_varmain, subTest_varc1,
                     subTest_varc2, sig_alpha, ttype):

    if np.linalg.norm(subTest_varc1) < 1e-6:
        diff_resid1 = subTest_diff
        numV = subTest_diff.shape[0]
    else:
        subTest_varc12 = sp.zeros((subTest_varc1.shape[0], 2))
        for i in range(subTest_varc1.shape[0]):
            subTest_varc12[i, 0] = subTest_varc1[i]
            subTest_varc12[i, 1] = subTest_varc2[i]
        print('regressing out 2 covariates')
        diff_resid1 = np.zeros(subTest_diff.shape)
        numV = subTest_diff.shape[0]
        for nv in tqdm(range(numV)):
            diff_resid1[nv, :] = multiLinReg_resid(subTest_varc12,
                                                   subTest_diff[nv, :])

    print('computing correlation against main variable')
    rval = np.zeros(numV)
    pval = np.zeros(numV)
    for nv in tqdm(range(numV)):
        if ttype == 'linear':
            rval[nv], pval[nv] = sp.stats.pearsonr(
                subTest_varmain, diff_resid1[nv, :])
        if ttype == 'group':
            rval[nv], pval[nv] = sp.stats.ttest_ind(
                subTest_varmain, diff_resid1[nv, :])
            #g = set(subTest_varmain)
            # for gnum in len(g):

    p = np.zeros(len(pval))
    p[pval <= sig_alpha] = 1
    p[np.isnan(pval)] = 0

    a = spio.loadmat('supp_data/USCBrain_grayordinate_labels.mat')
    labs = a['labels'].squeeze()
    labs[np.isnan(labs)] = 0
    labs[np.isnan(pval)] = 0
    pval_fdr = sp.zeros(numV)
    _, pv = fdrcorrection(pval[labs > 0], alpha=sig_alpha)
    pval_fdr[labs > 0] = pv

    pf = np.zeros(len(pval))
    pf[pval_fdr <= sig_alpha] = 1
    pf[labs == 0] = 0

    pval_fdr[labs == 0] = float("NAN")

    msg = str(np.sum(p)) + ' significant voxels found. ' + str(
        np.sum(pf)) + ' significant voxels found after FDR correction.'
    print(msg)

    return rval, pval, pval_fdr, msg


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
