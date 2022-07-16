%clear all;close all;clc;
addpath(genpath('/home/ajoshi/projects/bfp/src'))
fmri = '/home/ajoshi/sub08001/func/sub08001_rest_bold';
subbasename = '/home/ajoshi/sub08001/anat/sub08001_T1w';
config.FSLPATH = '/home/ajoshi/webware/fsl'
config.FSLOUTPUTTYPE='NIFTI_GZ';
config.AFNIPATH = '/home/ajoshi/abin';
config.FSLRigidReg=0;
config.MultiThreading=1;
BFPPATH='/home/ajoshi/projects/bfp';
[~,subid] = fileparts(subbasename);
sessionid = 'rest';

%function gen_alff_gord()
get_alff_gord(config, fmri, BFPPATH, subid, sessionid, subbasename);
