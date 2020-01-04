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
from statsmodels.stats.multitest import fdrcorrection

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
subIDs_grp1 = subIDs[group == 1]
sub_fname_grp1 = sub_fname[group == 1]

subIDs_grp2 = subIDs[group == 1]
sub_fname_grp2 = sub_fname[group == 1]

#%% makes file list for subcjects

tscore, pval = randpair_groupdiff(sub_fname_grp1,
                                  sub_fname_grp2,
                                  num_pairs=100,
                                  len_time=235)
##%%

vis_grayord_sigcorr(pval,
                    tscore,
                    cf.outname,
                    cf.out_dir,
                    int(cf.smooth_iter),
                    cf.save_surfaces,
                    cf.save_figures,
                    'True',
                    bfp_path=cf.bfp_path,
                    fsl_path=cf.fsl_path)

pval[sp.isnan(pval)] = .5

_, pval_fdr = fdrcorrection(pval)

vis_grayord_sigcorr(pval_fdr,
                    tscore,
                    cf.outname + '_fdr',
                    cf.out_dir,
                    int(cf.smooth_iter),
                    cf.save_surfaces,
                    cf.save_figures,
                    'False',
                    bfp_path=cf.bfp_path,
                    fsl_path=cf.fsl_path)

write_text_timestamp(log_fname,
                     'BFP Group difference pairwise analysis complete')
