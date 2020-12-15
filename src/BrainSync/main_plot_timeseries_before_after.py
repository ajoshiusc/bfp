# ||AUM||
import scipy.io as spio
import scipy as sp
import numpy as np
from fmri_methods_sipi import rot_sub_data
from surfproc import view_patch_vtk, patch_color_attrib
from dfsio import readdfs
import os
import glob
from brainsync import normalizeData, brainSync
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
import h5py
from tqdm import tqdm
from matplotlib.pyplot import figure,gcf
import matplotlib

# %%
# This is the name of the directory that contains the grayordinate files
#p_dir = '/big_disk/ajoshi/ADHD_Peking_gord'
#lst = glob.glob(p_dir + '/9*.filt.mat')
p_dir = '/data_disk/HCP_data/data'
lst = os.listdir(p_dir)

count1 = 0
# Get number of subjects
lst=lst[:6]
nsub = len(lst)
lst2 = list([])

for sub in tqdm(lst):
    fname = os.path.join(p_dir,sub,sub+'.rfMRI_REST1_LR.reduce3.ftdata.NLM_11N_hvar_25.mat')
    if not os.path.isfile(fname):
        continue
    df = spio.loadmat(fname)  # , 'r')
    d = df['ftdata_NLM'].T
    #   dataR = df['dataR']
    #   d = sp.concatenate((dataL, dataR), axis=1)
    d, _, _ = normalizeData(d)
    if count1 == 0:
        sub_data = np.zeros((d.shape[0], d.shape[1], nsub))

    sub_data[:, :, count1] = d
    lst2.append(fname)
    count1 += 1

# %% Compute pairwise distance
nSub = count1
sub_data = sub_data[:, :, :nSub]

print(nSub)
dist_all_orig = np.zeros([nSub, nSub])
dist_all_rot = dist_all_orig.copy()
#sub_data_orig = sub_data.copy()

ind1 =1
plt.rcParams.update({'font.size': 22})

plt.figure(num=1, figsize=(30, 3), dpi=300, facecolor='w', edgecolor='k')
plt.ion()
for ind2 in range(nSub):
    dist_all_orig[ind1, ind2] = sp.linalg.norm(sub_data[:, :, ind1] -
                                                sub_data[:, :, ind2])
    sub_data_rot, _ = brainSync(
        X=sub_data[:, :, ind1], Y=sub_data[:, :, ind2])
    dist_all_rot[ind1, ind2] = sp.linalg.norm(sub_data[:, :, ind1] -
                                                sub_data_rot)
    plt.plot(sub_data[:230,10000])
    print(ind1, ind2, dist_all_rot[ind1, ind2])
    plt.draw()
plt.savefig('before.png')

plt.show()
