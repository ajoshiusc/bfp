clc;clear all;close all;restoredefaultpath;
addpath(genpath('/project/ajoshi_27/code_farm/bfp/src'));
% Configuration file

%lst=dir('/home/rcf-40/ajoshi/aaj/Beijing_Zang/sub*');
lst= dir('/scratch1/ajoshi/parkinsons_fmri/neurocon/sub*');

for jj=1:length(lst)
    subbasename=['/scratch1/ajoshi/parkinsons_fmri/neurocon_bfp_out/',lst(jj).name,'/anat/',lst(jj).name];
    % % Outputs will be saved in the studydir/subid for each subject
    
    if ~exist([subbasename,'_T1w.SCT.GOrd.mat'],'file')
        %lst(jj).name
        fprintf(['sbatch --export=subid=',lst(jj).name,' /project/ajoshi_27/code_farm/bfp/supp_data/parkinson_slurm_neurocon.sh\n']);
    end
    
    
end


