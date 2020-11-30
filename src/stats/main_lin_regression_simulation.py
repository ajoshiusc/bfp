import scipy.io as spio
import scipy as sp
import numpy as np
import os
import sys

sys.path.append('../BrainSync')

from stats_utils import randpairs_regression_simulation
from surfproc import view_patch_vtk, patch_color_attrib, smooth_surf_function, smooth_patch
from dfsio import readdfs
from brainsync import normalizeData, brainSync
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import fdrcorrection
from stats_utils import randpairs_regression
from dev_utils import read_fcon1000_data
from grayord_utils import visdata_grayord, vis_grayord_sigpval
# ### Set the directorie
# s for the data and BFP software
from tqdm import tqdm
import time
# In[2]:

BFPPATH = '/ImagePTE1/ajoshi/code_farm/bfp'
FSL_PATH = '/usr/share/fsl/5.0'
# study directory where all the grayordinate files lie
DATA_DIR = '/big_disk/ajoshi/ADHD_Peking_gord'
CSV_FILE = '/data_disk/ADHD/ADHD_Peking_bfp/Peking_all_phenotypic.csv'

# ### Read CSV file to read the group IDs. This study has three subgroups:
# 1. Normal controls,
# 2. ADHD-hyperactive, and
# 3. ADHD-inattentive.

LEN_TIME = 235  # length of the time series
NUM_SUB = 150  # Number of subjects for the study


def main():

    print('Reading subjects')

    _, reg_var, sub_files = read_fcon1000_data(
        csv_fname=CSV_FILE,
        data_dir=DATA_DIR,
        reg_var_name='ADHD Index',  #'Verbal IQ',  #  #
        num_sub=NUM_SUB)

    # Shuffle reg_var and subjects for testing
    reg_var = sp.random.permutation(reg_var)
    #ran_perm = sp.random.permutation(len(reg_var))
    #reg_var = reg_var
    #sub_files = [sub_files[i] for i in range(len(reg_var))]

    t0 = time.time()
    print('performing stats based on random pairwise distances')

    corr_pval_max, corr_pval_fdr = randpairs_regression_simulation(
        bfp_path=BFPPATH,
        sub_files=sub_files,
        reg_var=reg_var,
        num_pairs=2000,  # 19900,
        nperm=2000,
        len_time=LEN_TIME,
        num_proc=4,
        pearson_fdr_test=False)
    t1 = time.time()

    print(t1 - t0)
    sp.savez(
        'pval_num_pairs2000_nsub200_nperm2000_simulation.npz',
        corr_pval_max=corr_pval_max,
        corr_pval_fdr=corr_pval_fdr)
    # corr_pval_max=a['corr_pval_max']
    # corr_pval_fdr=a['corr_pval_fdr']
    vis_grayord_sigpval(
        corr_pval_max,
        surf_name='rand_dist_corr_perm_pairs2000_max_simulation',
        out_dir='.',
        smooth_iter=1000,
        bfp_path=BFPPATH,
        fsl_path=FSL_PATH,
        sig_alpha=0.05)
    vis_grayord_sigpval(
        corr_pval_fdr,
        surf_name='rand_dist_corr_perm_pairs2000_fdr_simulation',
        out_dir='.',
        smooth_iter=1000,
        bfp_path=BFPPATH,
        fsl_path=FSL_PATH,
        sig_alpha=0.05)

    print('Results saved')


if __name__ == "__main__":
    main()