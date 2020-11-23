import time
from grayord_utils import save2volbord_bci
from dev_utils import read_fcon1000_data
from stats_utils import randpairs_regression
import scipy as sp
import numpy as np
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
        num_sub=NUM_SUB,
        gord=0)

    # Shuffle reg_var and subjects for testing
    #reg_var = sp.random.permutation(reg_var)
    #ran_perm = sp.random.permutation(len(reg_var))
    #reg_var = reg_var
    #sub_files = [sub_files[i] for i in range(len(reg_var))]

    t0 = time.time()
    print('performing stats based on random pairwise distances')

    corr_pval_max, corr_pval_fdr = randpairs_regression(
        bfp_path=BFPPATH,
        sub_files=sub_files,
        reg_var=reg_var,
        num_pairs=2000,  # 19900,
        nperm=2000,
        len_time=LEN_TIME,
        num_proc=6,
        pearson_fdr_test=False)
    t1 = time.time()

    print(t1 - t0)
    sp.savez(
        'pval_num_pairs2000_nsub200_nperm2000_bord.npz',
        corr_pval_max=corr_pval_max,
        corr_pval_fdr=corr_pval_fdr)
    # corr_pval_max=a['corr_pval_max']
    # corr_pval_fdr=a['corr_pval_fdr']

    save2volbord_bci((0.05-corr_pval_fdr)*np.float32(corr_pval_fdr < 0.05),
                     'pval_bord_ADHD_fdr_smooth1.5.nii.gz', bfp_path=BFPPATH, smooth_std=1.5)
    save2volbord_bci((0.05-corr_pval_max)*np.float32(corr_pval_max < 0.05),
                     'pval_bord_ADHD_max_smooth1.5.nii.gz', bfp_path=BFPPATH, smooth_std=1.5)
    save2volbord_bci((0.05-corr_pval_fdr)*np.float32(corr_pval_fdr < 0.05),
                     'pval_bord_ADHD_fdr_smooth1.nii.gz', bfp_path=BFPPATH, smooth_std=1.0)
    save2volbord_bci((0.05-corr_pval_max)*np.float32(corr_pval_max < 0.05),
                     'pval_bord_ADHD_max_smooth1.nii.gz', bfp_path=BFPPATH, smooth_std=1.0)
    save2volbord_bci((0.05-corr_pval_fdr)*np.float32(corr_pval_fdr < 0.05),
                     'pval_bord_ADHD_fdr_smooth0.5.nii.gz', bfp_path=BFPPATH, smooth_std=0.5)
    save2volbord_bci((0.05-corr_pval_max)*np.float32(corr_pval_max < 0.05),
                     'pval_bord_ADHD_max_smooth0.5.nii.gz', bfp_path=BFPPATH, smooth_std=0.5)

    print('Results saved')


if __name__ == "__main__":
    main()
