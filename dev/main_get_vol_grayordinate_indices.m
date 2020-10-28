% This file is a helper function. It takes grayordinates for the volume
% part from one sample cifti file. They grayordinates are specified in 2mm
% MNI space. They are transferred to the BCI space. For this purpose, first
% we coregistred BrainSuiteAtlas1 to BCI atlas and then used the svreg map.
clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/ImagePTE1/ajoshi/code_farm/bfp/dev/gifti-1.6'));
addpath(genpath('/ImagePTE1/ajoshi/code_farm/svreg/src'));
addpath(genpath('/ImagePTE1/ajoshi/code_farm/svreg/3rdParty'));
addpath(genpath('/big_disk/ajoshi/freesurfer/matlab'));
addpath(genpath('/ImagePTE1/ajoshi/code_farm/bfp/dev/Washington-University-cifti-matlab-27383b8'));
v=ft_read_cifti('/data_disk/HCP5/100307/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii');
bci=load_nii_BIG_Lab('/home/ajoshi/BrainSuite19b/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.bfc.nii.gz');
SZ=size(bci.img);

voxc=((v.transform)\([v.pos,ones(length(v.pos),1)]'))'; % This is in voxel coordinates
voxc=voxc(:,1:3);
save('../supp_data/MNI2mm_gord_vol_coord.mat','voxc')
voxc = (voxc)*2; % this is because the grayordinate voxel space was for 2mm res
voxc=voxc-1;
map=load_nii_BIG_Lab('/big_disk/ajoshi/bst1_processed/BrainSuiteAtlas1/mri.svreg.map.nii.gz');
mapx=map.img(:,:,:,1);
mapy=map.img(:,:,:,2);
mapz=map.img(:,:,:,3);
ind=~isnan(voxc(:,1));

% Mask of Grayordinates in BrainSuiteAtlas1 coordinates (MNI space)
v1=load_nii_BIG_Lab('/big_disk/ajoshi/bst1_processed/BrainSuiteAtlas1/mri.bfc.nii.gz');
id1=sub2ind(size(mapx),voxc(ind,1),voxc(ind,2),voxc(ind,3));
v1.hdr.dime.datatype=2;v1.img=0*v1.img;
v1.img(id1)=255;
save_untouch_nii_gz(v1,'../supp_data/Vol_MRI_grayord32k.mask.nii.gz');


% Generate mask of Grayordinates in BCI-DNI coordinates
mp1 = max(1,round(mapx(id1)));
mp2 = max(1,round(mapy(id1)));
mp3 = max(1,round(mapz(id1)));

bci_vol_ind=nan(length(voxc),1);
bci_vol_ind(ind)=sub2ind(SZ(1:3),mp1,mp2,mp3);

voxc(ind,:) = [mp1,mp2,mp3];
save('../supp_data/bci_gord_vol_coord.mat','voxc')
save('../supp_data/bci_grayordinates_vol_ind.mat','bci_vol_ind');

% Mask of Grayordinates in BCI-DNI coordinates (BCI space)
v=load_nii_BIG_Lab('/ImagePTE1/ajoshi/code_farm/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.nii.gz');
v.hdr.dime.datatype=2;v.img=0*v.img;
v.img(bci_vol_ind(ind))=255;
save_untouch_nii_gz(v,'../supp_data/Vol_BCI_grayord32k.mask.nii.gz');

