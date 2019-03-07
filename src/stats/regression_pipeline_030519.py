#%%
### Required Inputs
### Set the directories for BFP software
import sys
BFPPATH = '/home/sychoi/Documents/MATLAB/bfp/'
### Set directories for data
NDim = 5 # Dimensionality reduction for analysis
data_dir = '/NCAdisk/BFPtest/hypoxia_data' #input directory
out_dir = '/home/sychoi/Dropbox/SCD/Analysis/BOLD/022119/hypoxia_BFPstatTest' #output directory
csvfname = '/home/sychoi/Dropbox/SCD/Analysis/BOLD/022119/hypoxia/demohypoxia.csv' #csv file with demographics
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
from tqdm import tqdm
import os
import csv
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import fdrcorrection
import time

sys.path.append(BFPPATH + "/src/BrainSync/") 
from surfproc import view_patch_vtk, patch_color_attrib, smooth_surf_function, smooth_patch
from dfsio import readdfs,writedfs
from brainsync import IDrefsub_BrainSync, normalizeData, brainSync, groupBrainSync
sys.path.append(BFPPATH + "/src/stats/") 
from stats_utils import dist2atlas, load_bfp_data, read_demoCSV, sync2atlas, dist2atlas_reg, lin_reg, vis_save_pval, ref_avg_atlas, randpairsdist_reg_parallel, randpairsdist_reg


#%% 
print(out_dir + ": writing output directory")
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
#%% read demographic csv file
sub_ID, sub_fname, subAtlas_idx, reg_var, reg_cvar1, reg_cvar2 = read_demoCSV(csvfname,
                data_dir,
                file_ext,
                colsubj,
                colvar_exclude,
                colvar_atlas,
                colvar_main,
                colvar_reg1,
                colvar_reg2)

#%% ID reference subjects
print('Identifying subjects for atlas creation...')
subAtlas_fname = []; subAtlas_IDs = []
for ind in range(len(sub_ID)):
    if int(subAtlas_idx[ind]) !=0:
        sub = sub_ID[ind]
        #print(sub)
        fname = sub_fname[ind]
        subAtlas_fname.append(fname)  
        subAtlas_IDs.append(sub)
np.savetxt(out_dir + "/subjects_atlas.csv", subAtlas_IDs, delimiter=",", fmt='%s')
#%%
print('Identifying subjects for hypothesis testing...')
subTest_fname = []; subTest_IDs = []; 
for ind in range(len(sub_ID)):
    if int(subAtlas_idx[ind]) !=1:
        sub = sub_ID[ind]
        #print(sub)
        fname = sub_fname[ind]
        subTest_fname.append(fname)
        subTest_IDs.append(sub)  
np.savetxt(out_dir + "/subjects_testing.csv", subTest_IDs, delimiter=",", fmt='%s')
#%%
# reads reference data and creates atlas by Group BrainSync algorithm
dataAtlas = load_bfp_data(subAtlas_fname, LenTime)
refAtlas, q = IDrefsub_BrainSync(dataAtlas)
print(str(q))
#X2, Os, Costdif, TotalError = groupBrainSync(dataAtlas)

spio.savemat(os.path.join(out_dir + '/atlas.mat'), {'refAtlas': refAtlas})
del dataAtlas# Os, Costdif, TotalError

# reads subject data
dataTest = load_bfp_data(subTest_fname, LenTime)
#%%
# sync and calculates geodesic distances
syn_data = sync2atlas(refAtlas, dataTest)
#spio.savemat(os.path.join(out_dir + '/sync2atlas.mat'), {'syn_data': syn_data})
diff = dist2atlas(refAtlas, syn_data)
spio.savemat(os.path.join(out_dir + '/dist2atlas.mat'), {'diff': diff})
del dataTest
#%%
rcorr = sp.zeros(diff.shape[0])
r = sp.array(reg_var[:diff.shape[1]])
r = sp.absolute(r - sp.mean(r))
pcorr = sp.zeros(diff.shape[0])
for nv in range(diff.shape[0]):
    rho, pval  = sp.stats.pearsonr(diff[nv,:],r)
    rcorr[nv] = rho
    pcorr[nv] = pval
    
print(nv)
#%%
vis_save_pval(BFPPATH, pcorr, 'test', out_dir, smooth_iter=1500)

#%%
lsurf = readdfs(BFPPATH + '/supp_data/bci32kleft.dfs')
rsurf = readdfs(BFPPATH + '/supp_data/bci32kright.dfs')
a = spio.loadmat(BFPPATH + '/supp_data/USCBrain_grayord_labels.mat')
labs = a['labels']

lsurf.attributes = np.zeros((lsurf.vertices.shape[0]))
rsurf.attributes = np.zeros((rsurf.vertices.shape[0]))
lsurf = smooth_patch(lsurf, iterations=1500)
rsurf = smooth_patch(rsurf, iterations=1500)
labs[sp.isnan(labs)] = 0
print(pcorr.shape, labs.shape)
pcorr = pcorr*(labs > 0)

nVert = lsurf.vertices.shape[0]

p = abs(pcorr)

lsurf.attributes = p.squeeze()
lsurf.attributes = lsurf.attributes[:nVert]
rsurf.attributes = p.squeeze()
rsurf.attributes = rsurf.attributes[nVert:2*nVert]
lsurf = patch_color_attrib(lsurf, clim=[0, .05])
rsurf = patch_color_attrib(rsurf, clim=[0, .05])
lsurf.vColor[lsurf.attributes < 0, :] = .05
rsurf.vColor[rsurf.attributes < 0, :] = .05
writedfs(out_dir + '/rcorr.dfs',rsurf)
writedfs(out_dir + '/lcorr.dfs',lsurf)
spio.savemat(os.path.join(out_dir + '/rcorr.mat'), {'rcorr': rcorr})

print(lsurf.attributes.shape, nVert, lsurf.vColor.shape)
view_patch_vtk(lsurf, azimuth=100, elevation=180, roll=90, outfile=out_dir + '/l1corr.png', show=1)
view_patch_vtk(rsurf, azimuth=-100, elevation=180, roll=-90,
               outfile=out_dir + '/r1corr.png', show=1)
