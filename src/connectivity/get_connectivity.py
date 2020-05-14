import sys
import os
import numpy as np
sys.path.append('../stats')

#from surfproc import patch_color_attrib, smooth_surf_function
from brainsync import normalizeData

#from dfsio import readdfs, writedfs
from scipy import io as spio

#%%


def get_connectivity(data, labels, label_ids):
    #%%
    if type(data) == str:
        df = spio.loadmat(data)
        data = df['dtseries'].T

    num_time = data.shape[0]

    num_rois = len(label_ids)

    rtseries = np.zeros((num_time, num_rois))

    for i, id in enumerate(label_ids):

        idx = labels == id

        rtseries[:, i] = np.mean(data[:, idx], axis=1)

    rtseries, _, _ = normalizeData(rtseries)

    conn = np.corrcoef(rtseries.T)

    return conn


if __name__ == "__main__":

    BFPPATH = '/home/ajoshi/coding_ground/bfp'
    BrainSuitePath = '/home/ajoshi/BrainSuite19b/svreg'
    NDim = 31

    p_dir = '/ImagePTE1/ajoshi/fitbir/preproc/maryland_rao_v1/TBI_INVZV163RWK/BFP/TBI_INVZV163RWK/func/'
    sub = 'TBI_INVZV163RWK'
    atlas_labels = '/home/ajoshi/projects/bfp/supp_data/USCBrain_grayordinate_labels.mat'

    atlas = spio.loadmat(atlas_labels)

    gord_labels = atlas['labels'].squeeze()

    label_ids = np.unique(gord_labels)  # unique label ids

    # remove WM label from connectivity analysis
    label_ids = np.setdiff1d(label_ids, 2000)

    fname = os.path.join(p_dir, sub + '_rest_bold.32k.GOrd.mat')
    conn = get_connectivity(fname, gord_labels, label_ids)

    input('press any key')