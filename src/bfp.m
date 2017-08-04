%% BFP: BrainSuite fMRI Pipeline

%% User Inputs

t1 = '/home/ajoshi/coding_ground/bfp/data/sub-01_T1w.nii.gz';
fmri{1} = '/home/ajoshi/coding_ground/bfp/data/sub-01_ses-movie_task-movie_run-1_bold.nii.gz';
subbasename = '~/coding_ground/bfp/data/sub-01-run1/sub-01-run1';
BrainSuitePath='/home/ajoshi/BrainSuite17a';
bst_exe=fullfile('/home/ajoshi/coding_ground/bfp/src/cortical_extraction_nobse.sh');
svreg_exe=fullfile(BrainSuitePath,'svreg/svreg.sh');
BCIbasename=fullfile(BrainSuitePath,'svreg/BCI-DNI_brain_atlas/BCI-DNI_brain');

GOrdSurfIndFile='/home/ajoshi/coding_ground/bfp/dev/bci_grayordinates_surf_ind.mat';
GOrdVolIndFile='/home/ajoshi/coding_ground/bfp/dev/bci_grayordinates_vol_ind.mat';

RMFLAG=1;

fprintf('OS:%s\n',computer);
if ~strcmp(computer,'GLNXA64')
    error('OS is not supported!')
end

fprintf('Processing T1=%s\n',t1)
for ind=1:length(fmri)
    fprintf('Processing fMRI=%s\n',fmri{ind})
end

%% Create Directory Structure
fprintf('Creating Directory Structure\n');
subdir=fileparts(subbasename);
if RMFLAG && exist(subdir,'dir')
    rmdir(subdir,'s')
end

mkdir(subdir)
anatDir=fullfile(subdir,'anat');
fprintf('Creating Dir:%s\n',anatDir);
mkdir(anatDir);
t1new=fullfile(anatDir,'orig_hires.nii.gz');
t1ds=fullfile(anatDir,'orig.nii.gz'); 
copyfile(t1,t1new);
cmd=sprintf('flirt -in %s -ref %s -out %s -applyisoxfm 1',t1new,t1new,t1ds);
unix(cmd);
for ind = 1:length(fmri)
    outdir=fullfile(subdir,sprintf('func-%d',ind));
    fprintf('Creating Dir:%s\n',outdir);
    mkdir(outdir);    
    copyfile(fmri{ind},fullfile(outdir,'fmri.nii.gz'));
end
%% Skull Strip MRI

fprintf('Performing Skull Extraction\n');
bse=fullfile(BrainSuitePath,'bin','bse');
bseout=fullfile(anatDir,'orig.bse.nii.gz');
cmd=sprintf('%s --auto -i %s -o %s',bse,t1ds,bseout);
unix(cmd);
 
%% Coregister t1 to MNI Space
bsenew=fullfile(anatDir,'t1.nii.gz');
BSA=fullfile(BrainSuitePath,'svreg/BrainSuiteAtlas1/mri.bfc.nii.gz');
cmd=sprintf('flirt -ref %s -in %s -out %s',BSA,bseout,bsenew);
unix(cmd);
bsemask=fullfile(anatDir,'t1.mask.nii.gz');
cmd=sprintf('fslmaths %s -thr 0 -bin -mul 255 %s -odt char',bsenew,bsemask);
unix(cmd);
bsenew2=fullfile(anatDir,'t1.bse.nii.gz');
copyfile(bsenew,bsenew2);
 
%% Run BrainSuite
subbasename=fullfile(anatDir,'t1');
cmd=sprintf('%s %s',bst_exe,subbasename);
unix(cmd)
cmd=sprintf('%s %s %s',svreg_exe,subbasename,BCIbasename);
unix(cmd);

% 
%% Run Batch_Process Pipeline
% 
%% Transfer data to surface and then to USCBrain atlas surface and produce surface grayordinates
fprintf('Transferring data from subject to atlas\n');
for ind = 1:length(fmri)
    outdir=fullfile(subdir,sprintf('func-%d',ind));
    fmri2surfFile=fullfile(outdir,'fmri2surf.mat');
    GOrdsurfFile=fullfile(outdir,'fmri2surf_GOrd.mat');
    resample2surf(subbasename,fullfile(outdir,'fmri.nii.gz'),fmri2surfFile);
    generateSurfGOrdfMRI(GOrdSurfIndFile,fmri2surfFile,GOrdsurfFile);
end

%% Transfer data to volumetric grayordinates
fprintf('Generating Volumetric Grayordinates\n');
for ind = 1:length(fmri)
    outdir=fullfile(subdir,sprintf('func-%d',ind));
    GOrdVolFile=fullfile(outdir,'fmri2Vol_GOrd.mat');
    generateVolGOrdfMRI(GOrdVolIndFile,subbasename,fullfile(outdir,'fmri.nii.gz'),GOrdVolFile);
    combineSurfVolGOrdfMRI(GOrdSurfIndFile,GOrdVolFile);
end


%% Combine Surface and Volumetric Grayordinates and write gifti file
%