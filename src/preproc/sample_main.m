clc;clear all;close all;restoredefaultpath;
addpath(genpath('/home/rcf-40/ajoshi/aaj/git_sandbox/bfp'));
% Configuration file
configfile='/home/rcf-40/ajoshi/aaj/git_sandbox/bfp/supp_data/hpcconfig.ini';
%T1 and fMRI images in NIFTI-1 format (nii.gz)
 t1='/home/rcf-40/ajoshi/aaj/Beijing_Zang/sub01018/anat/mprage_anonymized.nii.gz';
 fmri='{''/home/rcf-40/ajoshi/aaj/Beijing_Zang/sub01018/func/rest.nii.gz''}';
% % Outputs will be saved in the studydir/subid for each subject
 studydir='/home/rcf-40/ajoshi/aaj/git_sandbox/bfp/data/';
 subid='sub01018';
 % fMRI files will be saved in BIDS format with sessionid in fmri names
 sessionid='{''rest''}';
 TR='2';
 
 bfpcommand = sprintf('disp %s %s %s %s %s %s %s', configfile, t1, fmri, studydir, subid, sessionid, TR);
cmd1=sprintf('/usr/usc/matlab/default/bin/matlab -nodisplay -nosplash -r "addpath(genpath(''/home/rcf-proj2/aaj/git_sandbox/svreg-matlab/src''));%s;exit;"',bfpcommand);

fprintf(cmd1);
unix(cmd1);

