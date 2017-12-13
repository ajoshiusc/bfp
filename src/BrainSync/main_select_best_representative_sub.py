# ||AUM||
import scipy.io
import scipy as sp
import numpy as np
from fmri_methods_sipi import rot_sub_data
from surfproc import view_patch_vtk, patch_color_attrib
from dfsio import readdfs
import os
from sklearn.manifold import MDS
import matplotlib.pyplot as plt

p_dir = '/big_disk/ajoshi/HCP_data/data'
p_dir_ref = '/big_disk/ajoshi/HCP_data'
lst = os.listdir(p_dir)

r_factor = 3
ref_dir = os.path.join(p_dir_ref, 'reference')
nClusters = 30

ref = '196750'#'100307'
print(ref + '.reduce' + str(r_factor) + '.LR_mask.mat')
fn1 = ref + '.reduce' + str(r_factor) + '.LR_mask.mat'
fname1 = os.path.join(ref_dir, fn1)
msk = scipy.io.loadmat(fname1)  # h5py.File(fname1);
dfs_right = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.\
a2009s.32k_fs.reduce3.right.dfs'))
dfs_right_sm = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.\
a2009s.32k_fs.reduce3.very_smooth.right.dfs'))
count1 = 0
rho_rho = []
rho_all = []
labs_all = sp.zeros((len(dfs_right.labels), len(lst)))

for sub in lst:
    data = scipy.io.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST1_LR.\
reduce3.ftdata.NLM_11N_hvar_25.mat'))
    LR_flag = msk['LR_flag']
    LR_flag = np.squeeze(LR_flag) == 0
    data = data['ftdata_NLM']
    temp = data[LR_flag, :]
    m = np.mean(temp, 1)
    temp = temp - m[:, None]
    s = np.std(temp, 1)+1e-16
    temp = temp/s[:, None]
    d = temp
    if count1 == 0:
        sub_data = sp.zeros((d.shape[0], d.shape[1], len(lst)))

    sub_data[:, :, count1] = d
    count1 += 1
    print count1,

nSub = sub_data.shape[2]
rperm = sp.random.permutation(dfs_right_sm.vertices.shape[0])
#rperm=range(dfs_right_sm.vertices.shape[0])
dist_all_orig = sp.zeros([nSub, nSub])
dist_all_rot = dist_all_orig.copy()
#sub_data[:,:,1]=sub_data[rperm,:,1]
sub_data_orig = sub_data.copy()

for ind1 in range(nSub):
    for ind2 in range(nSub):
        dist_all_orig[ind1, ind2] = sp.linalg.norm(sub_data_orig[:, :, ind1] -
                                                   sub_data_orig[:, :, ind2])
        sub_data_rot, _ = rot_sub_data(ref=sub_data[:, :, ind1],
                                       sub=sub_data[:, :, ind2])
        dist_all_rot[ind1, ind2] = sp.linalg.norm(sub_data[:, :, ind1] -
                                                  sub_data_rot)
        print ind1, ind2


sp.savez('rot_pairwise_dist_all_sub_by_sub.npz', dist_all_rot=dist_all_rot,
         dist_all_orig=dist_all_orig, lst=lst)
######

a = sp.load('rot_pairwise_dist_all_sub_by_sub.npz')
q = sp.argmin(a['dist_all_rot'].sum(1))
m=MDS(n_components=3,dissimilarity='precomputed')
e=m.fit_transform(a['dist_all_rot'])
print(e)
fig, ax = plt.subplots()
ax.scatter(e[:,0],e[:,1])
for i in range(e.shape[0]):
    ax.annotate(lst[i], (e[i,0],e[i,1]))
