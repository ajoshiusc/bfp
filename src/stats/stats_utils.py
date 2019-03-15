""" This module contains helpful utility function for running statistics using BFP """

import csv
import os
import scipy as sp
import numpy as np
import scipy.io as spio
from tqdm import tqdm
import itertools
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
from sklearn.decomposition import PCA
import sklearn
from surfproc import view_patch_vtk, patch_color_attrib, smooth_surf_function, smooth_patch
from dfsio import readdfs, writedfs
import sys
sys.path.append('../BrainSync')
import multiprocessing
from functools import partial
from brainsync import normalizeData, brainSync


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

def corr_pearson_fdr(X, Y, nperm=1000):
    #X: nsub x vertices
    #Y: cognitive scores nsub X 1
    num_vert = X.shape[1]

    corr_pval = sp.zeros(num_vert)
    for ind in tqdm(range(num_vert)):
        _, corr_pval[ind] = sp.stats.pearsonr(X[:,ind], Y)

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