#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 19:50:23 2019

@author: sychoi
"""

import os
import sys
import scipy.io as spio
import numpy as np

bfp_path = '/home/sychoi/Documents/MATLAB/bfp/'
sys.path.append(os.path.join(bfp_path, 'src/stats/') )
from read_data_utils import load_bfp_data, write_text_timestamp
sys.path.append(os.path.join(str(bfp_path), 'src/BrainSync/')) 
from brainsync import IDrefsub_BrainSync, groupBrainSync, brainSync
#%%
dirname = "/NCAdisk/SCD_structural_analysis/BOLD_study/SCD_BOLDdata/"
out_dir = '/NCAdisk/SCD_structural_analysis/BOLD_study/BOLD_Analysis/032519/hyperoxia/'
ext = '_hyperoxia_bold.32k.GOrd.filt.mat'
LenTime = 150
#%% checks and loads data
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
flog = os.path.join(out_dir, 'bfp_stat_log.txt')
write_text_timestamp(flog, 'starting group-MDS_atlas')

os.chdir(bfp_path)
subj = os.listdir(dirname)
sub_fname = []
sub_ID = []
n=0

for i in range(len(subj)):
    print(str(n))
    n=n+1
    fname = os.path.join(dirname, subj[i], 'func', subj[i] + ext)
    if os.path.isfile(fname):
        df = spio.loadmat(fname)
        data = df['dtseries'].T
        if int(data.shape[0]) == LenTime:
            sub_ID.append(subj[i])
            sub_fname.append(fname)
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