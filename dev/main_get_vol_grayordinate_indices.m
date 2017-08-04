% This file is a helper function. It takes grayordinates for the volume
% part from one sample cifti file. They grayordinates are specified in 2mm
% MNI space. They are transferred to the BCI space. For this purpose, first
% we coregistred BrainSuiteAtlas1 to BCI atlas and then used the svreg map.
clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/bfp/dev/gifti-1.6'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/3rdParty'));
addpath(genpath('/big_disk/ajoshi/freesurfer/matlab'));
addpath(genpath('/home/ajoshi/coding_ground/bfp/dev/Washington-University-cifti-matlab-27383b8'));
v=ft_read_cifti('/big_disk/ajoshi/HCP5/100307/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii');
bci=load_nii_BIG_Lab('/home/ajoshi/BrainSuite17a/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.bfc.nii.gz');
SZ=size(bci.img);

voxc=((v.transform)\([v.pos,ones(length(v.pos),1)]'))'; % This is in voxel coordinates
voxc = voxc*2; % this is because the grayordinate voxel space was for 2mm res

map=load_nii_BIG_Lab('/home/ajoshi/Desktop/BrainSuiteAtlas1/mri.svreg.map.nii.gz');
mapx=map.img(:,:,:,1);
mapy=map.img(:,:,:,2);
mapz=map.img(:,:,:,3);
ind=~isnan(voxc(:,1));
voxc(ind,1) = max(min(voxc(ind,1),SZ(1)),1);
voxc(ind,2) = max(min(voxc(ind,2),SZ(2)),1);
voxc(ind,3) = max(min(voxc(ind,3),SZ(3)),1);

id1=sub2ind(size(mapx),voxc(ind,1),voxc(ind,2),voxc(ind,3));
mp1 = max(1,round(mapx(id1)));
mp2 = max(1,round(mapy(id1)));
mp3 = max(1,round(mapz(id1)));

bci_vol_ind=nan(length(voxc),1);
bci_vol_ind(ind)=sub2ind(SZ(1:3),mp1,mp2,mp3);

save('bci_grayordinates_vol_ind.mat','bci_vol_ind');



