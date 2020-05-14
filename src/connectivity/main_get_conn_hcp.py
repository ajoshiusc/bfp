import glob
import os
import sys

import h5py
import numpy as np
#from dfsio import readdfs, writedfs
from scipy import io as spio
from tqdm import tqdm

#from surfproc import patch_color_attrib, smooth_surf_function
sys.path.append('../stats')

from brainsync import normalizeData
from get_connectivity import get_connectivity

BFPPATH = '/home/ajoshi/coding_ground/bfp'
BrainSuitePath = '/home/ajoshi/BrainSuite19b/svreg'
NDim = 31

#%%

if __name__ == "__main__":

    atlas_labels = '/home/ajoshi/projects/bfp/supp_data/USCBrain_grayordinate_labels.mat'

    atlas = spio.loadmat(atlas_labels)

    gord_labels = atlas['labels'].squeeze()

    label_ids = np.unique(gord_labels)  # unique label ids

    # remove WM label from connectivity analysis
    label_ids = np.setdiff1d(label_ids, (2000, 0))

    p_dir = '/data_disk/HCP_100_Andrew_Preprocessed_05_13'
    lst = glob.glob(p_dir + '/*')
    count1 = 0

    # Get number of subjects
    nsub = len(lst)

    conn_mat = np.zeros((len(label_ids), len(label_ids), nsub))

    for subno, subdir in tqdm(enumerate(lst)):
        fname = os.path.join(subdir, 'rfMRI_1_LR.mat')
        if not os.path.isfile(fname):
            continue

        f = h5py.File(fname, 'r')


        # Convert Andrew's format to grayordinate format
        dataL = np.array(f['dataL'])
        dataR = np.array(f['dataR'])
        dataS = np.array(f['dataS'])
        d = np.concatenate((dataL, dataR, dataS), axis=1)

        # Get connectivity matrix for each subject
        conn_mat[:, :, subno] = get_connectivity(d,
                                                 labels=gord_labels,
                                                 label_ids=label_ids)

    input('press any key')
