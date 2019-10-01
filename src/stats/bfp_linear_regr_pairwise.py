#%%
config_file = '/home/ajoshi/coding_ground/bfp/src/stats/sample_config_stats.ini'
#%%#%%
### Import the required librariesimport configparser
import sys
import os
from stats_utils import randpairs_regression, multiLinReg_resid
from grayord_utils import vis_grayord_sigpval
import scipy.io as spio
import scipy as sp
import numpy as np
import configparser
from tqdm import tqdm

### Import BrainSync libraries
config = configparser.ConfigParser()
config.read(config_file)
section = config.sections()
bfp_path = config.get('inputs', 'bfp_path')
from read_data_utils import load_bfp_data, read_demoCSV, write_text_timestamp, readConfig
os.chdir(bfp_path)
cf = readConfig(config_file)
from brainsync import IDrefsub_BrainSync, groupBrainSync, generate_avgAtlas
from stats_utils import dist2atlas, sync2atlas, multiLinReg_corr
from grayord_utils import vis_grayord_sigcorr

#%%
log_fname = os.path.join(cf.out_dir, 'bfp_linregr_stat_log.txt')
write_text_timestamp(log_fname, 'Config file used: ' + config_file)
if not os.path.isdir(cf.out_dir):
    os.makedirs(cf.out_dir)
write_text_timestamp(log_fname,
                     "All outputs will be written in: " + cf.out_dir)
# read demographic csv file
sub_ID, sub_fname, subAtlas_idx, reg_var, reg_cvar1, reg_cvar2 = read_demoCSV(
    cf.csv_fname, cf.data_dir, cf.file_ext, cf.colsubj, cf.colvar_exclude,
    cf.colvar_atlas, cf.colvar_main, cf.colvar_reg1, cf.colvar_reg2)
#%% makes file list for subjects
print('Identifying subjects for atlas creation and hypothesis testing...')
subTest_fname = []
subTest_IDs = []
subAtlas_fname = []
subAtlas_IDs = []
for ind in range(len(sub_ID)):
    sub = sub_ID[ind]
    fname = sub_fname[ind]
    if int(subAtlas_idx[ind]) == 1:
        subAtlas_fname.append(fname)
        subAtlas_IDs.append(sub)
    else:
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
    if cf.test_all == 'False':
        if int(subAtlas_idx[ind]) != 1:
            count1 += 1
    else:
        count1 += 1
    subTest_varmain[count1] = varmain
    subTest_varc1[count1] = varc1
    subTest_varc2[count1] = varc2

del sub_ID, sub_fname, subAtlas_idx, reg_var, reg_cvar1, reg_cvar2, fname, sub, count1, ind, numT

import csv
with open(cf.out_dir + "/subjects_testing.csv", 'w') as csvfile:
    csv.writer(csvfile).writerows(
        zip(subTest_IDs, subTest_varmain, subTest_varc1, subTest_varc2))
np.savetxt(
    cf.out_dir + "/subjects_atlas.csv", subAtlas_IDs, delimiter=",", fmt='%s')
write_text_timestamp(
    log_fname,
    str(len(subAtlas_IDs)) + ' subjects will be used for atlas creation.')
write_text_timestamp(
    log_fname,
    str(len(subTest_IDs)) + ' subjects will be used for hypothesis testing.')

#%%
# reads reference data and creates atlas by BrainSync algorithm
if len(cf.atlas_fname) != 0:
    write_text_timestamp(
        log_fname,
        'User Option: User defined atlas will be used ' + cf.atlas_fname)
    df = spio.loadmat(cf.atlas_fname)
    atlas_data = df['atlas_data']
    del df
else:
    subAtlas_data = load_bfp_data(subAtlas_fname, int(cf.lentime))
    if cf.atlas_groupsync == 'True':
        write_text_timestamp(
            log_fname,
            'User Option: Group BrainSync algorithm will be used for atlas creation'
        )
        atlas_data, _, _, _ = groupBrainSync(subAtlas_data)
        write_text_timestamp(log_fname, 'Done creating atlas')
    else:
        write_text_timestamp(
            log_fname,
            'User Option: representative subject will be used for atlas creation'
        )
        subRef_data, subRef_num = IDrefsub_BrainSync(subAtlas_data)
        write_text_timestamp(
            log_fname, 'Subject number ' + str(subAtlas_IDs[subRef_num]) +
            ' will be used for atlas creation')
        atlas_data = generate_avgAtlas(subRef_data, subAtlas_data)
    del subAtlas_data
    spio.savemat(
        os.path.join(cf.out_dir + '/atlas.mat'), {'atlas_data': atlas_data})

#%% Compute pairwise distance and perform regression
corr_pval_max, corr_pval_fdr = randpairs_regression(
    bfp_path=cf.bfp_path,
    sub_files=subTest_fname,
    reg_var=subTest_varmain,
    num_pairs=5000,  # 19900,
    nperm=2000,
    len_time=int(cf.lentime),
    num_proc=6,
    pearson_fdr_test=False)

# subTest_data = load_bfp_data(subTest_fname, int(cf.lentime))
# subTest_syndata = sync2atlas(atlas_data, subTest_data)
# subTest_diff = dist2atlas(atlas_data, subTest_syndata)
# spio.savemat(
#     os.path.join(cf.out_dir + '/dist2atlas.mat'),
#     {'subTest_diff': subTest_diff})
# del subTest_data, subTest_syndata
#%% computes correlation after controlling for two covariates
""" rval, pval, pval_fdr, msg = multiLinReg_corr(subTest_diff, subTest_varmain,
                                             subTest_varc1, subTest_varc2)
write_text_timestamp(log_fname, 'Done runnning linear regression. ' + msg)
 """#%%

vis_grayord_sigpval(
    corr_pval_max,
    surf_name=cf.outname + '_max',
    out_dir=cf.out_dir,
    smooth_iter=int(cf.smooth_iter),
    bfp_path=cf.bfp_path,
    fsl_path=cf.fsl_path)

vis_grayord_sigpval(
    corr_pval_fdr,
    surf_name=cf.outname + 'fdr',
    out_dir=cf.out_dir,
    smooth_iter=int(cf.smooth_iter),
    bfp_path=cf.bfp_path,
    fsl_path=cf.fsl_path)

write_text_timestamp(log_fname, 'BFP regression analysis complete')