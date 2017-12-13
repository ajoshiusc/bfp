# ||AUM||
import scipy.io
import scipy as sp
import numpy as np
from fmri_methods_sipi import rot_sub_data
from surfproc import view_patch_vtk, patch_color_attrib
from dfsio import readdfs
import os
from brainsync import normalizeData, brainSync
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
import h5py

p_dir = '/big_disk/ajoshi/HCP_Data_For_Haleh'
lst = os.listdir(p_dir)
count1 = 0
fn1 = 'rfMRI_1_LR.mat'
fname1 = os.path.join(p_dir, fn1)
mskf = h5py.File(fname1, 'r')
print(list(mskf.keys()))
mskf.close()
nsub = 0
for sub in lst:
    fname = os.path.join(p_dir, sub, 'rfMRI_1_LR.mat')
    if not os.path.isfile(fname):
        continue
    nsub += 1

for sub in lst:
    fname = os.path.join(p_dir, sub, 'rfMRI_1_LR.mat')
    if not os.path.isfile(fname):
        continue
    df = h5py.File(fname, 'r')
    dataL = df['dataL']
    dataR = df['dataR']
    d = sp.concatenate((dataL, dataR), axis=1)
    d, _, _ = normalizeData(d)
    if count1 == 0:
        sub_data = sp.zeros((d.shape[0], d.shape[1], nsub))

    sub_data[:, :, count1] = d
    count1 += 1
    print count1,


nSub = sub_data.shape[2]
dist_all_orig = sp.zeros([nSub, nSub])
dist_all_rot = dist_all_orig.copy()
sub_data_orig = sub_data.copy()

for ind1 in range(nSub):
    for ind2 in range(nSub):
        dist_all_orig[ind1, ind2] = sp.linalg.norm(sub_data_orig[:, :, ind1] -
                                                   sub_data_orig[:, :, ind2])
        sub_data_rot, _ = brainSync(X=sub_data[:, :, ind1],
                                    Y=sub_data[:, :, ind2])
        dist_all_rot[ind1, ind2] = sp.linalg.norm(sub_data[:, :, ind1] -
                                                  sub_data_rot)
        print ind1, ind2


sp.savez('Haleh_pairwise_dist_all_sub_by_sub.npz', dist_all_rot=dist_all_rot,
         dist_all_orig=dist_all_orig, lst=lst)
######

a = sp.load('rot_pairwise_dist_all_sub_by_sub.npz')
q = sp.argmin(a['dist_all_rot'].sum(1))
m = MDS(n_components=3, dissimilarity='precomputed')
e = m.fit_transform(a['dist_all_rot'])
print(e)
fig, ax = plt.subplots()
ax.scatter(e[:, 0], e[:, 1])
for i in range(e.shape[0]):
    ax.annotate(lst[i], (e[i, 0], e[i, 1]))
