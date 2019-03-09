#%%
### Required Inputs
useGroupSync = False # False if you'd like to create reference atlas by identifying one representative subject
### Set the directories for BFP software
import sys
BFPPATH = '/home/sychoi/Documents/MATLAB/bfp/'
### Set directories for data
NDim = 5 # Dimensionality reduction for analysis
data_dir = '/NCAdisk/BFPtest/hypoxia_data' #input directory
out_dir = '/home/sychoi/Dropbox/SCD/Analysis/BOLD/022119/hypoxia_BFPstatTest' #output directory
csvfname = '/home/sychoi/Dropbox/SCD/Analysis/BOLD/022119/hypoxia/demohypoxia.csv' #csv file with demographics
outname = 'hgbCorr-agesexregress' # file subnames for result outputs (example: outdir/outname_pval.png)
colsubj = 'SubjectID'
colvar_main = 'HGB'
colvar_reg1 = 'age'
colvar_reg2 = 'Sex'
colvar_exclude = 'Exclude'
colvar_atlas = 'Reference'
file_ext = '_hypoxia_bold.32k.GOrd.filt.mat' #input file extension
LenTime = 150 #number of timepoints
#%%
### Import the required libraries
import scipy.io as spio
import scipy as sp
import numpy as np
import os
### Import BrainSync libraries
sys.path.append(BFPPATH + "/src/BrainSync/") 
from brainsync import IDrefsub_BrainSync, groupBrainSync, generate_avgAtlas
sys.path.append(BFPPATH + "/src/stats/") 
from stats_utils import LinReg_corr, dist2atlas, load_bfp_data, read_demoCSV, sync2atlas, vis_save_pval

#%% 
os.chdir(BFPPATH)
print(out_dir + ": writing output directory")
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
# read demographic csv file
sub_ID, sub_fname, subAtlas_idx, reg_var, reg_cvar1, reg_cvar2 = read_demoCSV(csvfname,
                data_dir,
                file_ext,
                colsubj,
                colvar_exclude,
                colvar_atlas,
                colvar_main,
                colvar_reg1,
                colvar_reg2)
#%% makes file list for subjects
print('Identifying subjects for atlas creation and hypothesis testing...')
subTest_fname = []; subTest_IDs = []; subAtlas_fname = []; subAtlas_IDs = []
for ind in range(len(sub_ID)):
    sub = sub_ID[ind]
    fname = sub_fname[ind]
    if int(subAtlas_idx[ind]) !=1:        
        subTest_fname.append(fname)
        subTest_IDs.append(sub)
    else:
        subAtlas_fname.append(fname)  
        subAtlas_IDs.append(sub)
len(subTest_IDs)

count1=0
subTest_varmain = sp.zeros(len(subTest_IDs)); subTest_varc1 = sp.zeros(len(subTest_IDs)); subTest_varc2 = sp.zeros(len(subTest_IDs))
for ind in range(len(sub_ID)):
    varmain = reg_var[ind]
    varc1 = reg_cvar1[ind]
    varc2 = reg_cvar2[ind]
    if int(subAtlas_idx[ind]) !=1:  
        subTest_varmain[count1] = varmain
        subTest_varc1[count1] = varc1
        subTest_varc2[count1] = varc2
        count1 +=1

np.savetxt(out_dir + "/subjects_testing.csv", subTest_IDs, delimiter=",", fmt='%s')
np.savetxt(out_dir + "/subjects_atlas.csv", subAtlas_IDs, delimiter=",", fmt='%s')
print(str(len(subAtlas_IDs)) + ' subjects will be used for atlas creation.')
print(str(len(subTest_IDs)) + ' subjects will be used for hypothesis testing.')
#%%
# reads reference data and creates atlas by BrainSync algorithm
subAtlas_data = load_bfp_data(subAtlas_fname, LenTime)
if useGroupSync == True:
    print('User Option: Group BrainSync algorithm will be used for atlas creation')
    atlas_data, _, _, _ = groupBrainSync(subAtlas_data)
else:
    print('User Option: representative subject will be used for atlas creation')
    subRef_data, subRef_num = IDrefsub_BrainSync(subAtlas_data)
    print('Subject number ' + str(subRef_num) + ' will be used for atlas creation')
    atlas_data = generate_avgAtlas(subRef_data, subAtlas_data)
    
spio.savemat(os.path.join(out_dir + '/atlas.mat'), {'atlas_data': atlas_data})
del subAtlas_data
#%% sync and calculates geodesic distances
subTest_data = load_bfp_data(subTest_fname, LenTime)
subTest_syndata = sync2atlas(atlas_data, subTest_data)
subTest_diff = dist2atlas(atlas_data, subTest_syndata)
spio.savemat(os.path.join(out_dir + '/dist2atlas.mat'), {'subTest_diff': subTest_diff})
del subTest_data, subTest_syndata
#%% computes correlation after controlling for two covariates
rval, pval, pval_fdr = LinReg_corr(subTest_diff, subTest_varmain, subTest_varc1, subTest_varc2 )
#%%
vis_save_pval(BFPPATH, pval, outname , out_dir, smooth_iter=1000)
vis_save_pval(BFPPATH, pval_fdr, outname + '_fdr', out_dir, smooth_iter=1000)