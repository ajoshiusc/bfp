""" This module contains helpful utility function for running statistics using BFP """

import csv
import os
import scipy as sp
import scipy.io as spio
from brainsync import normalizeData, brainSync
from tqdm import tqdm
from statsmodels.stats.multitest import fdrcorrection
from sklearn.decomposition import PCA
import statsmodels.api as sm


def read_fcon1000_data(csv_fname,
                       data_dir,
                       reg_var_name='Verbal IQ',
                       num_sub=5,
                       len_time=250):
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

            # Load data and normalize it
            data = spio.loadmat(fname)
            data = data['dtseries'].T
            data, _, _ = normalizeData(data)

            if count1 == 0:
                sub_data = sp.zeros((len_time, data.shape[1], num_sub))

            # Truncate the data at a given number of time samples This is needed because
            # BrainSync needs same number of time sampples
            sub_data[:, :, count1] = data[:len_time, ]
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

    return sub_ids, reg_var, sub_data


def sync2atlas(atlas, sub_data):
    print('Syncing to atlas, assume that the data is normalized')

    # Assume that the sub_data is already normalized
    syn_data = sp.zeros(sub_data.shape)
    for ind in tqdm(range(sub_data.shape[2])):
        syn_data[:, :, ind], _ = brainSync(X=atlas, Y=sub_data[:, :, ind])

    return syn_data


def dist2atlas_reg(ref_atlas, sub_data, reg_var):
    """ Perform regression stats based on distance to atlas """
    print('dist2atlas_reg, assume that the data is normalized')

    num_vert = sub_data.shape[1]
    num_sub = sub_data.shape[2]

    diff = sp.zeros([sub_data.shape[1], num_sub])

    # Compute distance to atlas
    for ind in tqdm(range(num_sub)):
        Y2, _ = brainSync(X=ref_atlas, Y=sub_data[:, :, ind])
        diff[:, ind] = sp.sum((Y2 - ref_atlas)**2, axis=0)

    corr_pval = sp.zeros(num_vert)
    for v in tqdm(range()):
        _, corr_pval[v] = sp.stats.pearsonr(diff[v, :], reg_var)

    corr_pval[sp.isnan(corr_pval)] = .5
    _, corr_pval_fdr = fdrcorrection(corr_pval)

    return corr_pval, corr_pval_fdr


def lin_reg(ref_atlas, sub_data, reg_var, ndim=20):
    """ Perform regression stats based on distance to atlas """

    num_vert = sub_data.shape[1]
    num_sub = sub_data.shape[2]

    print('Computing PCA basis function from the atlas')
    pca = PCA(n_components=ndim)
    pca.fit(ref_atlas.T)

    rData = sp.zeros((ndim, sub_data.shape[1], num_sub))
    for ind in tqdm(range(num_sub)):
        Y2, _ = brainSync(X=ref_atlas, Y=sub_data[:, :, ind])
        rData[:, :, ind] = pca.transform(Y2.T).T

    pval_linreg = sp.zeros(num_vert)

    for v in tqdm(range(num_vert)):
        X = rData[:, v, :]
        X = sm.add_constant(X.T)
        est = sm.OLS(reg_var, X)
        pval_linreg[v] = est.fit().f_pvalue

    print('Regression is done')

    pval_linreg[sp.isnan(pval_linreg)] = .5
    _, pval_linreg_fdr = fdrcorrection(pval_linreg)

    return pval_linreg, pval_linreg_fdr

