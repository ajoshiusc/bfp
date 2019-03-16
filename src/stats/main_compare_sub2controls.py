# coding: utf-8

# # Regression Pipeline

# This notebook serves as a template for performing a group difference comparison using BFP and BrainSync. The steps in this pipeline can be easily customized to suite your study. Here, we use data from ADHD200 dataset available through http://fcon_1000.projects.nitrc.org/indi/adhd200/. Specifically, we use the Peking dataset. We will correlate cognitive score and fmri signal deviation from normals. This will find the areas that are associated with fMRI signal.
#
# The pipeline is written in Python (Jupyter Notebook). We assume thet BrainSuite and BFP are installed on your computer. Install the required python libraries listed below in the script. We recommend using Anaconda python distribution.
#
# The steps for running the group comparison are:
#
# * Process the fMRI and T1 data of subjects using BFP.
# * Set the paths in group analysis script.
# * Run the regression script.
#
# As an input, we assume that all the subjects data has been preprocessed using BFP. Specifically, we will use the grayordinate data produced by BFP. Also a CSV file containing group labels is assume as an input.
#
# First, we use a set of normal control subjects to build an average atlas. Currently, it is done by using BrainSync to synchronize all subject data to an individual and then averaging the synchronized data. In the future, we can use group brainsync included in BFP.
#
# For a population of subjects, we will compute correlation between norm of synchronized fMRI signal of subjects and atlas, and the cognitive scores of the subjects.
#
# In the future, we will compute multivariate regression.
#
# ### Import the required libraries

# In[1]:

import scipy.io as spio
import scipy as sp
import numpy as np
from surfproc import view_patch_vtk, patch_color_attrib, smooth_surf_function, smooth_patch
from dfsio import readdfs
import os
import sys
sys.path.append('../BrainSync')
from brainsync import normalizeData, brainSync
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import fdrcorrection
from stats_utils import compare_sub2ctrl, read_gord_data, dist2atlas_reg, lin_reg, vis_save_pval, randpairsdist_reg_parallel, randpairsdist_reg
# ### Set the directories for the data and BFP software
from tqdm import tqdm
import time
# In[2]:

BFPPATH = '/home/ajoshi/coding_ground/bfp'

# study directory where all the grayordinate files lie
CTRL_DIR = '/big_disk/ajoshi/CN_new_gord'

SUB_DATA = '/big_disk/ajoshi/for_cleveland/bfpout/study13072/func/study13072_rest_bold.32k.GOrd.filt.mat'

LEN_TIME = 235  # length of the time series
NUM_CTRL = 100  # Number of control subjects for the study


def main():

    print('Reading subjects')

    ctrl_files = read_gord_data(data_dir=CTRL_DIR, num_sub=NUM_CTRL)

    t0 = time.time()
    print('performing stats based on random pairwise distances')

    pval = compare_sub2ctrl(
        bfp_path=BFPPATH,
        sub_file=SUB_DATA,
        ctrl_files=ctrl_files,
        num_pairs=2000,
        nperm=1000,
        len_time=LEN_TIME,
        num_proc=4,
        fdr_test=False)
    t1 = time.time()

    print(t1 - t0)

    vis_save_pval(
        bfp_path=BFPPATH,
        out_dir='.',
        pval_map=pval,
        surf_name='rand_dist_corr_par')

    print('Results saved')


if __name__ == "__main__":
    main()