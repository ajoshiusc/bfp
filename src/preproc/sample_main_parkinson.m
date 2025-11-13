clc;clear all;close all;restoredefaultpath;
%addpath(genpath('/big_disk/ajoshi/coding_ground/bfp/supp_data'))
addpath(genpath('/project/ajoshi_27/code_farm/bfp/src'));
%    1050345 rest 2

% Set the input arguments
configfile='/project/ajoshi_27/code_farm/bfp/supp_data/config_bfp_preproc_minimal_carc.ini';


sublist = dir('/scratch1/ajoshi/parkinsons_fmri/taowu/sub*');

for j=1:length(sublist)

    subid= sublist(j).name;

    t1=['/scratch1/ajoshi/parkinsons_fmri/taowu/',subid,'/anat/',subid,'_T1w.nii.gz'];
    fmri=['/scratch1/ajoshi/parkinsons_fmri/taowu/',subid,'/func/',subid,'_task-resting_bold.nii.gz'];

    %fmri='/home/ajoshi/Downloads/BFP_issues/ACTL005/ACTL005.BOLD.resting.nii.gz';
    studydir='/scratch1/ajoshi/parkinsons_fmri/taowu_bfp_out';
    sessionid='task-resting';
    TR='';

    setenv('BrainSuiteMCR','/project/ajoshi_27/MATLAB_Runtime/v912');
    try
        % Call the bfp function
        bfp(configfile, t1, fmri, studydir, subid, sessionid,TR);
    catch
        fprintf('***** sub with error: %s\n ***** \n',subid);
    end

end
