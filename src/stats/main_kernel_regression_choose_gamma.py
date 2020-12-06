import time
from dev_utils import read_fcon1000_data
from ker_reg_utils import kernel_regression_choose_gamma
import numpy as np
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
NUM_SUB = 15  # Number of subjects for the study


def main():

    print('Reading subjects')

    _, reg_var, sub_files = read_fcon1000_data(
        csv_fname=CSV_FILE,
        data_dir=DATA_DIR,
        reg_var_name='ADHD Index',  # 'Verbal IQ',  #  #
        num_sub=100)

    # Shuffle reg_var and subjects for testing
    #reg_var = np.random.permutation(reg_var)
    #ran_perm = sp.random.permutation(len(reg_var))
    #reg_var = reg_var
    #sub_files = [sub_files[i] for i in range(len(reg_var))]

    t0 = time.time()
    print('performing stats based on kernel regression')
    reg_var = reg_var[15:15+NUM_SUB]
    sub_files = sub_files[15:15+NUM_SUB]

    gamma_max = kernel_regression_choose_gamma(
        bfp_path=BFPPATH,
        sub_files=sub_files,
        reg_var=reg_var,
        nperm=2000,
        len_time=LEN_TIME,
        num_proc=6,
        fdr_test=False)

    t1 = time.time()

    print(t1 - t0)

    print('optimum gamma is %s' % gamma_max)


if __name__ == "__main__":
    main()
