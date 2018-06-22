clc;clear all;close all;restoredefaultpath;
addpath(genpath('/home/rcf-proj2/aaj/git_sandbox/bfp'));
% Configuration file
configfile='/home/rcf-proj2/aaj/git_sandbox/bfp/supp_data/hpcconfig.ini';
%T1 and fMRI images in NIFTI-1 format (nii.gz)

lst=dir('/home/rcf-proj2/aaj/ADHD_Peking/data/Peking_all/*');

for jj=length(lst):-1:1
    if ~exist(['/home/rcf-proj2/aaj/ADHD_Peking_bfp/',lst(jj).name,'/func/',lst(jj).name,'_rest_bold.32k.GOrd.mat'],'file')
        unix(['sbatch --export=subid=',lst(jj).name,' /home/rcf-proj2/aaj/git_sandbox/bfp/supp_data/fcon1000_slurm.sh']);
    end
end

