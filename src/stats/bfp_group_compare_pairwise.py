#%%
config_file = '/home/ajoshi/coding_ground/bfp/src/stats/sample_config_stats_ADHD.ini'
#%%#%%
### Import the required librariesimport configparser
import sys
import os
import scipy.io as spio
import scipy as sp
import numpy as np
import configparser

### get BFP directory from config file
config = configparser.ConfigParser()
config.read(config_file)
section = config.sections()
bfp_path = config.get('inputs', 'bfp_path')
sys.path.append(os.path.join(bfp_path, 'src/stats/'))
from read_data_utils import load_bfp_data, read_demoCSV, write_text_timestamp, readConfig, read_demoCSV_list
os.chdir(bfp_path)
cf = readConfig(config_file)

### Import BrainSync libraries
from stats_utils import randpair_groupdiff
from grayord_utils import vis_grayord_sigcorr

#%%
log_fname = os.path.join(cf.out_dir, 'bfp_group_stat_log.txt')
write_text_timestamp(log_fname, 'Config file used: ' + config_file)
if not os.path.isdir(cf.out_dir):
    os.makedirs(cf.out_dir)
write_text_timestamp(log_fname,
                     "All outputs will be written in: " + cf.out_dir)
# read demographic csv file
subIDs, sub_fname, group, reg_var, reg_cvar1, reg_cvar2 = read_demoCSV(
    cf.csv_fname, cf.data_dir, cf.file_ext, cf.colsubj, cf.colvar_exclude,
    cf.colvar_group, cf.colvar_main, cf.colvar_reg1, cf.colvar_reg2)

group = np.int16(group)
# for boolan indexing, need to convert to numpy array
subIDs = np.array(subIDs)
sub_fname = np.array(sub_fname)

print('Identifying subjects for each group...')
subIDs_grp1 = subIDs[group == 0]
sub_fname_grp1 = sub_fname[group == 0]

subIDs_grp2 = subIDs[group == 1]
sub_fname_grp2 = sub_fname[group == 1]

#%% makes file list for subjects

randpair_groupdiff(sub_fname_grp1,
                   sub_fname_grp2,
                   num_pairs=10000,
                   len_time=235)
#
#pairs_grp1, pair_dist_grp1 = rand_pair_dist(sub_files=sub_fname_grp1,
#                                            num_pairs=num_pairs)

write_text_timestamp(
    log_fname,
    str(len(sub_ID)) + ' subjects will be used for hypothesis testing.')
#%%
subTest_data = load_bfp_data(subTest_fname, int(cf.lentime))
subTest_syndata = sync2atlas(atlas_data, subTest_data)
subTest_diff, _ = dist2atlas(atlas_data, subTest_syndata)
spio.savemat(os.path.join(cf.out_dir + '/dist2atlas.mat'),
             {'subTest_diff': subTest_diff})
del subTest_data, subTest_syndata

#%% computes correlation after controlling for two covariates
rval, pval, pval_fdr = multiLinReg_corr(subTest_diff, subTest_varmain,
                                        subTest_varc1, subTest_varc2)
#%%
vis_grayord_sigcorr(pval, rval, cf.outname, cf.out_dir, int(cf.smooth_iter),
                    cf.save_surfaces, cf.save_figures, 'True')
vis_grayord_sigcorr(pval, rval, cf.outname + '_fdr', cf.out_dir,
                    int(cf.smooth_iter), cf.save_surfaces, cf.save_figures,
                    'False')
write_text_timestamp(log_fname, 'BFP regression analysis complete')