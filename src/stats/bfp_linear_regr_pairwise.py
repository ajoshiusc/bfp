#%%
config_file = '/home/ajoshi/coding_ground/bfp/src/stats/sample_config_stats_ADHD.ini'
#%%#%%
### Import the required librariesimport configparser
import sys
import os
from stats_utils import randpairs_regression, multiLinReg_resid, LinReg_resid
from grayord_utils import vis_grayord_sigpval
import scipy.io as spio
import scipy as sp
import numpy as np
import configparser
from tqdm import tqdm
from sklearn.linear_model import LinearRegression

### Import BrainSync libraries
config = configparser.ConfigParser()
config.read(config_file)
section = config.sections()
bfp_path = config.get('inputs', 'bfp_path')
from read_data_utils import load_bfp_data, read_demoCSV, write_text_timestamp, readConfig
os.chdir(bfp_path)
cf = readConfig(config_file)
from stats_utils import multiLinReg_corr
from grayord_utils import vis_grayord_sigcorr

#%%
log_fname = os.path.join(cf.out_dir, 'bfp_linregr_stat_log.txt')
write_text_timestamp(log_fname, 'Config file used: ' + config_file)
if not os.path.isdir(cf.out_dir):
    os.makedirs(cf.out_dir)
write_text_timestamp(log_fname,
                     "All outputs will be written in: " + cf.out_dir)
# read demographic csv file
sub_ID, sub_fname, _, reg_var, reg_cvar1, reg_cvar2 = read_demoCSV(
    cf.csv_fname, cf.data_dir, cf.file_ext, cf.colsubj, cf.colvar_exclude,
    cf.colvar_atlas, cf.colvar_main, cf.colvar_reg1, cf.colvar_reg2)
#%% makes file list for subjects
print(
    'Identifying subjects for hypothesis testing, no atlas needs to be created...'
)
subTest_fname = []
subTest_IDs = []
for ind in range(len(sub_ID)):
    sub = sub_ID[ind]
    fname = sub_fname[ind]
    subTest_fname.append(fname)
    subTest_IDs.append(sub)

if cf.test_all == 'False':
    numT = len(subTest_IDs)
    write_text_timestamp(
        log_fname,
        "User Option: Only subjects not used for atlas creation will be used for hypothesis testing"
    )
else:
    numT = len(sub_ID)
    subTest_fname = sub_fname
    subTest_IDs = sub_ID
    write_text_timestamp(
        log_fname,
        "User Option: All subjects will be used for hypothesis testing")

count1 = -1
subTest_varmain = sp.zeros(numT)
subTest_varc1 = sp.zeros(numT)
subTest_varc2 = sp.zeros(numT)
for ind in range(len(sub_ID)):
    varmain = reg_var[ind]
    varc1 = reg_cvar1[ind]
    varc2 = reg_cvar2[ind]
    count1 += 1
    subTest_varmain[count1] = varmain
    subTest_varc1[count1] = varc1
    subTest_varc2[count1] = varc2

del sub_ID, sub_fname, reg_var, reg_cvar1, reg_cvar2, fname, sub, count1, ind, numT

import csv
with open(cf.out_dir + "/subjects_testing.csv", 'w') as csvfile:
    csv.writer(csvfile).writerows(
        zip(subTest_IDs, subTest_varmain, subTest_varc1, subTest_varc2))
write_text_timestamp(
    log_fname,
    ' There is no atlas needed for pairwise tests, all the subjects will be used for hypothesis testing.'
)
write_text_timestamp(
    log_fname,
    str(len(subTest_IDs)) + ' subjects will be used for hypothesis testing.')

#%% Do Linear regression on the covariates
write_text_timestamp(
    log_fname,
    'Doing linear regression on the covariates from the main variable.')

subTest_varc12 = sp.zeros((subTest_varc1.shape[0], 2))
subTest_varc12[:, 0] = subTest_varc1
subTest_varc12[:, 1] = subTest_varc2
regr = LinearRegression()
regr.fit(subTest_varc12, subTest_varmain)
#print(regr.coef_)
pre = regr.predict(subTest_varc12)
subTest_varmain2 = subTest_varmain - pre

#%% Compute pairwise distance and perform regression
corr_pval_max, corr_pval_fdr = randpairs_regression(
    bfp_path=cf.bfp_path,
    sub_files=subTest_fname,
    reg_var=subTest_varmain2,
    num_pairs=2000,  # 19900,
    nperm=2000,
    len_time=int(cf.lentime),
    num_proc=1,
    pearson_fdr_test=False)

#%% Visualization of the results
vis_grayord_sigpval(corr_pval_max,
                    surf_name=cf.outname + '_max',
                    out_dir=cf.out_dir,
                    smooth_iter=int(cf.smooth_iter),
                    bfp_path=cf.bfp_path,
                    fsl_path=cf.fsl_path)

vis_grayord_sigpval(corr_pval_fdr,
                    surf_name=cf.outname + 'fdr',
                    out_dir=cf.out_dir,
                    smooth_iter=int(cf.smooth_iter),
                    bfp_path=cf.bfp_path,
                    fsl_path=cf.fsl_path)

write_text_timestamp(log_fname, 'BFP regression analysis complete')