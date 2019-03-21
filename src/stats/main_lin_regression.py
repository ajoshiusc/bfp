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
from stats_utils import dist2atlas_reg, randpairsdist_reg_parallel, randpairsdist_reg
from dev_utils import read_fcon1000_data
from grayord_utils import visdata_grayord
# ### Set the directorie
# s for the data and BFP software
from tqdm import tqdm
import time
# In[2]:

BFPPATH = '/home/ajoshi/coding_ground/bfp'
FSL_PATH = '/usr/share/fsl/5.0'
# study directory where all the grayordinate files lie
DATA_DIR = '/deneb_disk/ADHD_Peking_gord'
CSV_FILE = '/deneb_disk/ADHD_Peking_bfp/Peking_all_phenotypic.csv'

# ### Read CSV file to read the group IDs. This study has three subgroups:
# 1. Normal controls,
# 2. ADHD-hyperactive, and
# 3. ADHD-inattentive.

LEN_TIME = 235  # length of the time series
NUM_SUB = 20  # Number of subjects for the study


def main():

    print('Reading subjects')

    _, reg_var, sub_files = read_fcon1000_data(
        csv_fname=CSV_FILE,
        data_dir=DATA_DIR,
        reg_var_name='ADHD Index',  #'Verbal IQ',  #  #
        num_sub=NUM_SUB)

    # Shuffle reg_var and subjects for testing
    # reg_var = sp.random.permutation(reg_var)
    #ran_perm = sp.random.permutation(len(reg_var))
    #    reg_var = reg_var[100:200]
    #    sub_files = [sub_files[i] for i in range(100, 200)]

    t0 = time.time()
    print('performing stats based on random pairwise distances')

    corr_pval = randpairsdist_reg_parallel(
        bfp_path=BFPPATH,
        sub_files=sub_files,
        reg_var=reg_var,
        num_pairs=5000,
        nperm=1000,
        len_time=LEN_TIME,
        num_proc=4,
        fdr_test=False)
    t1 = time.time()

    print(t1 - t0)

    visdata_grayord(
        corr_pval,
        surf_name='rand_dist_corr_vol',
        out_dir='.',
        smooth_iter=1000,
        colorbar_lim=[0, 1],
        colormap='jet_r',
        save_dfs=True,
        save_png=True,
        bfp_path=BFPPATH,
        fsl_path=FSL_PATH)

    print('Results saved')


if __name__ == "__main__":
    main()