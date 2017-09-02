clc;clear all;close all;restoredefaultpath;
addpath(genpath('/home/rcf-proj2/aaj/git_sandbox/bfp'));
% Configuration file
configfile='/home/rcf-proj2/aaj/git_sandbox/bfp/supp_data/hpcconfig.ini';
%T1 and fMRI images in NIFTI-1 format (nii.gz)

lst=dir('/home/rcf-proj2/aaj/Beijing_Zang/sub*');

for jj=1:length(lst)
   unix(['qsub -v subid=',lst(jj).name,' /home/rcf-proj2/aaj/git_sandbox/bfp/supp_data/fcon1000_qsub.sh &']);
end

