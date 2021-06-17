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
colin27=load_nii_BIG_Lab('/ImagePTE1/ajoshi/mni_colin27_1998/mni_colin27_BST/colin27_t1_tal_lin.bse.nii.gz');
SZ=size(colin27.img);

voxc=((v.transform)\([v.pos,ones(length(v.pos),1)]'))'; % This is in voxel coordinates
voxc=voxc(:,1:3);
save('../supp_data/colin27_2mm_gord_vol_coord.mat','voxc')
voxc = (voxc)*2; % this is because the grayordinate voxel space was for 2mm res
voxc=voxc-1;
% Mask of Grayordinates in Colin27 coordinates (MNI space)
v1=load_nii_BIG_Lab('/ImagePTE1/ajoshi/mni_colin27_1998/mni_colin27_BST/colin27_t1_tal_lin.bfc.nii.gz');
ind=~isnan(voxc(:,1));
id1=sub2ind(size(v1.img),voxc(ind,1),voxc(ind,2),voxc(ind,3));
v1.hdr.dime.datatype=2;v1.img=0*v1.img;
v1.img(id1)=255;
save_untouch_nii_gz(v1,'../supp_data/Vol_colin27_grayord32k.mask.nii.gz');


col_vol_ind=nan(length(voxc),1);
col_vol_ind(ind)=sub2ind(SZ(1:3),voxc(ind,1),voxc(ind,2),voxc(ind,3));

save('../supp_data/colin27_grayordinates_vol_ind.mat','col_vol_ind');

