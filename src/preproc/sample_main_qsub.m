clc;clear all;close all;restoredefaultpath;
addpath(genpath('/home/rcf-40/ajoshi/aaj/git_sandbox/bfp'));
% Configuration file
configfile='/home/rcf-40/ajoshi/aaj/git_sandbox/bfp/supp_data/hpcconfig.ini';
%T1 and fMRI images in NIFTI-1 format (nii.gz)

lst=dir('/home/rcf-40/ajoshi/aaj/Beijing_Zang/sub*');

for jj=1:length(lst)
    t1=['/home/rcf-40/ajoshi/aaj/Beijing_Zang/',lst(jj).name,'/anat/mprage_anonymized.nii.gz'];
    fmri=['{''/home/rcf-40/ajoshi/aaj/Beijing_Zang/',lst(jj).name,'/func/rest.nii.gz''}'];
    studydir='/home/rcf-40/ajoshi/aaj/git_sandbox/bfp/data/';
    subid=lst(jj).name;
    sessionid='{''rest''}';
    TR='2';
    
    bfpcommand = sprintf('bfp %s %s %s %s %s %s %s', configfile, t1, fmri, studydir, subid, sessionid, TR);
    cmd1=sprintf('/usr/usc/matlab/default/bin/matlab -nodisplay -nosplash -r "addpath(genpath(''/home/rcf-proj2/aaj/git_sandbox/svreg-matlab/src''));%s;exit;"',bfpcommand);
    
    fprintf(cmd1);
    unix(['qsub ',cmd1]);
    
end