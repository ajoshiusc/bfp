#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 19:50:23 2019

@author: sychoi
"""
#%%
config_file = '/home/sychoi/Dropbox/Vanderbilt/Analysis/HCP/stats/Analysis/2024-08-26/atlas_WMmask/initWMAtlas.ini'
#%%#%%
### Import the required librariesimport configparser
import sys
import os
import scipy.io as spio
import scipy as sp
import numpy as np
import configparser
import csv
from sklearn.linear_model import LinearRegression

### Import BrainSync libraries
config = configparser.ConfigParser()
config.read(config_file)
section = config.sections()
bfp_path = config.get('inputs','bfp_path')
sys.path.append(os.path.join(bfp_path, 'src/stats/') )
sys.path.append(os.path.join(str(bfp_path), 'src/BrainSync/'))
from read_data_utils import load_bfp_dataT, read_demoCSV, write_text_timestamp,readConfig,load_bfp_dataT_dist2atlas
os.chdir(bfp_path)
cf = readConfig(config_file)
from brainsync import IDrefsub_BrainSync, groupBrainSync, generate_avgAtlas
from stats_utils import dist2atlas, sync2atlas, multiLinReg_corr
from stats_utils import randpairs_regression, multiLinReg_resid, LinReg_resid, multiLinReg_corr
from grayord_utils import vis_grayord_sigcorr, vis_grayord_sigpval
#%%
if not os.path.isdir(cf.out_dir):
    os.makedirs(cf.out_dir)
log_fname = os.path.join(cf.out_dir, 'bfp_stat.log')
write_text_timestamp(log_fname, 'Config file used: ' + config_file +"\n All outputs will be written in: " + cf.out_dir )
# read demographic csv file
sub_ID, sub_fname, subAtlas_idx, reg_var, reg_cvar1, reg_cvar2 = read_demoCSV(cf.csv_fname,
                cf.data_dir,
                cf.file_ext,
                cf.colsubj,
                cf.colvar_exclude,
                cf.colvar_atlas,
                cf.colvar_main,
                cf.colvar_reg1,
                cf.colvar_reg2,
                cf.matcht,
                cf.lentime)
#%% makes file list for subjects

os.chdir(bfp_path)
subj = os.listdir(dirname)
sub_fname = []
sub_ID = []
n=0

pbar = tqdm(total=len(subj))
for i in range(len(subj)):
    print(str(n))
    n=n+1
    #fname = os.path.join(dirname, subj[i], 'func', subj[i] + ext)
    fname = os.path.join(dirname, subj[i], 'func', ext)
    if os.path.isfile(fname):
        df = spio.loadmat(fname)
        data = df['dtseries'].T
        if int(data.shape[0]) == LenTime:
            sub_ID.append(subj[i])
            sub_fname.append(fname)
    pbar.update(1)
np.savetxt(out_dir + "/subjects.csv", sub_ID, delimiter=",", fmt='%s')
write_text_timestamp(flog, 'data for ' + str(len(sub_ID)) + ' subjects will be used for atlas creation')
sub_data = load_bfp_data(sub_fname, LenTime)
#%% run group sync
groupatlas_data, _, _, _ = groupBrainSync(sub_data)
spio.savemat(os.path.join(out_dir + '/group_atlas.mat'), {'groupatlas_data': groupatlas_data})
write_text_timestamp(flog, 'completed group sync')
#%% find representative subject
subRef_data, q = IDrefsub_BrainSync(sub_data)
write_text_timestamp(flog, 'representative subject is ' + str(sub_ID[q]))
#%% sync group atlas to representative subject
atlas_data, _ = brainSync(subRef_data, groupatlas_data)
spio.savemat(os.path.join(out_dir + '/atlas.mat'), {'atlas_data': atlas_data})
write_text_timestamp(flog, 'group atlas synced to representative subject. group-MDS_atlas done.')