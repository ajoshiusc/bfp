clc;clear all;close all;restoredefaultpath;
%addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/supp_data'))
addpath(genpath('/home/ajoshi/Projects/bfp/src'));
%    1050345 rest 2

% Set the input arguments
configfile='/home/ajoshi/Projects/bfp/supp_data/config_bfp_preproc_minimal.ini';


sublist = dir('/deneb_disk/parkinsons_fmri/taowu/sub*');

for j=1:length(sublist)

    subid= sublist(j).name;

    t1=['/deneb_disk/parkinsons_fmri/taowu/',subid,'/anat/',subid,'_T1w.nii.gz'];
    fmri=['/deneb_disk/parkinsons_fmri/taowu/',subid,'/func/',subid,'_task-resting_bold.nii.gz'];

    %fmri='/home/ajoshi/Downloads/BFP_issues/ACTL005/ACTL005.BOLD.resting.nii.gz';
    studydir='/deneb_disk/parkinsons_fmri/taowu_bfp_out';
    sessionid='task-resting';
    TR='';

    setenv('BrainSuiteMCR','/home/ajoshi/Software/MATLAB/MATLAB_Runtime/R2023a');
    try
        % Call the bfp function
        bfp(configfile, t1, fmri, studydir, subid, sessionid,TR);
    catch
        fprintf('***** sub with error: %s\n ***** \n',subid);
    end

end
