import sys
import os
import numpy as np
sys.path.append('../stats')

#from surfproc import patch_color_attrib, smooth_surf_function
from brainsync import normalizeData

#from dfsio import readdfs, writedfs
from scipy import io as spio

BFPPATH = '/home/ajoshi/coding_ground/bfp'
BrainSuitePath = '/home/ajoshi/BrainSuite19b/svreg'
NDim = 31

#%%
p_dir = '/ImagePTE1/ajoshi/fitbir/preproc/maryland_rao_v1/TBI_INVZV163RWK/BFP/TBI_INVZV163RWK/func/'
sub = 'TBI_INVZV163RWK'

#%%
fname = os.path.join(p_dir, sub + '_rest_bold.32k.GOrd.mat')
df = spio.loadmat(fname)
data = df['dtseries'].T

num_time = data.shape[0]

d, _, _ = normalizeData(data)

atlas_labels = '/home/ajoshi/projects/bfp/supp_data/USCBrain_grayordinate_labels.mat'


atlas = spio.loadmat(atlas_labels)

labels = atlas['labels'].squeeze()
label_ids = np.unique(np.mod(labels,1000))

num_rois=len(label_ids)
Conn = np.zeros((num_rois,num_rois))

rtseries = np.zeros((num_time,num_rois))

for i, id in enumerate(label_ids):
    print(id)

    idx = labels==id

    rtseries[:,i] = np.mean(data[:,idx], axis=1)

rtseries, _, _ = normalizeData(rtseries)

conn = np.corrcoef(rtseries.T)


return conn, 
input("Press Enter to continue...")




