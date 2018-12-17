clc;clear all;close all;restoredefaultpath;
addpath(genpath('/home/rcf-40/ajoshi/aaj/git_sandbox/bfp'));
% Configuration file
configfile='/home/rcf-40/ajoshi/aaj/git_sandbox/bfp/supp_data/hpcconfig.ini';
%T1 and fMRI images in NIFTI-1 format (nii.gz)

lst=dir('/home/rcf-40/ajoshi/aaj/Beijing_Zang/sub*');

for jj=1:length(lst)
    subbasename=['/home/rcf-40/ajoshi/aaj/Beijing_Zang_bfp/',lst(jj).name,'/anat/',lst(jj).name];
    % % Outputs will be saved in the studydir/subid for each subject
    
    if ~exist([subbasename,'_T1w.SCT.GOrd.mat'],'file')
        lst(jj).name
        unix(['sbatch --export=subid=',lst(jj).name,' /home/rcf-proj2/aaj/git_sandbox/bfp/supp_data/fcon1000_slurm_orig.sh']);
    end
    
    
end

