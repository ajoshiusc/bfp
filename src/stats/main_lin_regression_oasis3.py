import scipy.io as spio
import scipy as sp
import numpy as np
from surfproc import view_patch_vtk, patch_color_attrib, smooth_surf_function, smooth_patch
from dfsio import readdfs
import os
import sys
sys.path.append('../BrainSync')
from brainsync import normalizeData, brainSync
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import fdrcorrection
from stats_utils import randpairs_regression
#from dev_utils import read_fcon1000_data
from utils_oasis3 import read_oasis3_data
from grayord_utils import visdata_grayord, vis_grayord_sigpval
# ### Set the directorie
# s for the data and BFP software
from tqdm import tqdm
import time

BFPPATH = '/home/ajoshi/projects/bfp'
FSL_PATH = '/usr/share/fsl/5.0'
# study directory where all the grayordinate files lie
DATA_DIR = '/deneb_disk/bfp_oasis3'
CSV_FILE = '/home/ajoshi/projects/AD_processing/src/oasis3_bfp.csv'

# ### Read CSV file to read the group IDs. This study has three subgroups:
# 1. Normal controls,
# 2. ADHD-hyperactive, and
# 3. ADHD-inattentive.

LEN_TIME = 164  # length of the time series
NUM_SUB = 350  # Number of subjects for the study


def main():

    print('Reading subjects')

    _, reg_var, sub_files = read_oasis3_data(
        csv_fname=CSV_FILE,
        data_dir=DATA_DIR,
        reg_var_name='UDSB9',  #'Verbal IQ',  #  #
        num_sub=NUM_SUB,
        len_time=LEN_TIME)

    # Shuffle reg_var and subjects for testing
    #reg_var = sp.random.permutation(reg_var)
    #ran_perm = sp.random.permutation(len(reg_var))
    #reg_var = reg_var
    #sub_files = [sub_files[i] for i in range(len(reg_var))]

    sub_files=sub_files #[50:100]
    reg_var = reg_var #[50:100]
    t0 = time.time()
    print('performing stats based on random pairwise distances')

    corr_pval_max, corr_pval_fdr = randpairs_regression(
        bfp_path=BFPPATH,
        sub_files=sub_files,
        reg_var=reg_var,
        num_pairs=20000,  # 19900,
        nperm=2000,
        len_time=LEN_TIME,
        num_proc=6,
        pearson_fdr_test=False)
    t1 = time.time()

    print(t1 - t0)
    sp.savez(
        'pval_num_pairs20000_nsub350_nperm2000.npz',
        corr_pval_max=corr_pval_max,
        corr_pval_fdr=corr_pval_fdr)
    # corr_pval_max=a['corr_pval_max']
    # corr_pval_fdr=a['corr_pval_fdr']
    vis_grayord_sigpval(
        corr_pval_max,
        surf_name='rand_dist_corr_perm_pairs20000_max',
        out_dir='.',
        smooth_iter=1000,
        bfp_path=BFPPATH,
        fsl_path=FSL_PATH,
        sig_alpha=0.05)
    vis_grayord_sigpval(
        corr_pval_fdr,
        surf_name='rand_dist_corr_perm_pairs20000_fdr',
        out_dir='.',
        smooth_iter=1000,
        bfp_path=BFPPATH,
        fsl_path=FSL_PATH,
        sig_alpha=0.05)

    print('Results saved')


if __name__ == "__main__":
    main()