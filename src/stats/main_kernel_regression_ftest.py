import time
from tqdm import tqdm
from grayord_utils import visdata_grayord, vis_grayord_sigpval
from dev_utils import read_fcon1000_data
from stats_utils import kernel_regression_ftest
from statsmodels.stats.multitest import fdrcorrection
from sklearn.decomposition import PCA
from brainsync import normalizeData, brainSync
import scipy.io as spio
import scipy as sp
import numpy as np
from surfproc import view_patch_vtk, patch_color_attrib, smooth_surf_function, smooth_patch
from dfsio import readdfs
import os
import sys
sys.path.append('../BrainSync')
# ### Set the directorie
# s for the data and BFP software
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
NUM_SUB = 50  # Number of subjects for the study


def main():

    print('Reading subjects')

    _, reg_var, sub_files = read_fcon1000_data(
        csv_fname=CSV_FILE,
        data_dir=DATA_DIR,
        reg_var_name='ADHD Index',  # 'Verbal IQ',  #  #
        num_sub=NUM_SUB)

    # Shuffle reg_var and subjects for testing
    reg_var = np.random.permutation(reg_var)
    #ran_perm = sp.random.permutation(len(reg_var))
    #reg_var = reg_var
    #sub_files = [sub_files[i] for i in range(len(reg_var))]
    #sub_files=sub_files[30:30+15]
    #reg_var = reg_var[30:30+15]

    t0 = time.time()
    print('performing stats based on kernel regression')
    
    pval_kr_ftest, pval_kr_ftest_fdr = kernel_regression_ftest(
        bfp_path=BFPPATH,
        sub_files=sub_files,
        reg_var=reg_var,
        nperm=2000,
        len_time=LEN_TIME,
        num_proc=4,
        fdr_test=False,
        simulation=False)
    t1 = time.time()

    print(t1 - t0)
    np.savez(
        'pval_KR_ftest.npz',
        pval=pval_kr_ftest,
        pval_fdr=pval_kr_ftest_fdr)
    # corr_pval_max=a['corr_pval_max']
    # corr_pval_fdr=a['corr_pval_fdr']
    
    #a=np.load('pval_KR_ftest5null.npz')
    #pval_kr_ftest = a['pval']
    vis_grayord_sigpval(
        pval_kr_ftest,
        surf_name='KR_pval_ftest',
        out_dir='.',
        smooth_iter=1000,
        bfp_path=BFPPATH,
        fsl_path=FSL_PATH,
        sig_alpha=0.5)
    vis_grayord_sigpval(
        pval_kr_ftest_fdr,
        surf_name='KR_pval_ftest_fdr',
        out_dir='.',
        smooth_iter=1000,
        bfp_path=BFPPATH,
        fsl_path=FSL_PATH,
        sig_alpha=0.05)

    print('Results saved')


if __name__ == "__main__":
    main()
