clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/ImagePTE1/ajoshi/code_farm/bfp/dev/gifti-1.6'));
addpath(genpath('/ImagePTE1/ajoshi/code_farm/svreg/src'));
addpath(genpath('/ImagePTE1/ajoshi/code_farm/svreg/3rdParty'));

BSTDir = '/home/ajoshi/BrainSuite19b';

usclobes_file = '/ImagePTE1/ajoshi/code_farm/svreg/USCLobes/BCI-DNI_brain.label.nii.gz';
target = fullfile(BSTDir, 'svreg','BCI-DNI_brain_atlas','BCI-DNI_brain.pvc.frac.nii.gz');

resampled_target = [tempname(),'.nii.gz'];

resampled_usclobes_file =  [tempname(),'.nii.gz'];
svreg_resample(usclobes_file, resampled_usclobes_file,'-res','3','3','3','nearest')
svreg_resample(target, resampled_target,'-res','3','3','3')


vv=load_nii_BIG_Lab(resampled_target);
ind = find(vv.img(:)>0);

usclobes = load_nii_BIG_Lab(resampled_usclobes_file);
usclobes_bord = usclobes.img(ind);

save('usclobes_bord.mat', 'usclobes_bord')
delete(resampled_target);
