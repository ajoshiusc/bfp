
% The purpose of this file is to generate a file that gives a vector of
% grayordinate labels, its id lookup tables
clc;clear all;close all;
addpath(genpath('/home/ajoshi/coding_ground/bfp/dev/gifti-1.6'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/3rdParty'));

% Read surface indices
load('../supp_data/bci_grayordinates_surf_ind.mat')
load('../supp_data/HCP_32k_Label.mat')
hcp_isnan=isnan(brainstructure(1:2*32492));

l= readdfs('/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.left.mid.cortex.dfs');

r= readdfs('/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.right.mid.cortex.dfs');

surf_labs=[l.labels(ind_left);r.labels(ind_right)];

surf_labs(hcp_isnan)=0;

% Read indices for volume and read labels from the corresponding volume

load('/home/ajoshi/coding_ground/bfp/supp_data/bci_grayordinates_vol_ind.mat')
bci_vol_ind(isnan(bci_vol_ind))=1;
vl = load_nii('/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.label.nii.gz');

vol_labs=vl.img(bci_vol_ind);

% Concatenate surface and volume labels
labels=[surf_labs;vol_labs(1+length(surf_labs):end)];

% Read XML and create lookup table
xml1 = xml2struct('/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/brainsuite_labeldescription.xml');

num_rois = length(xml1.labelset.label);
Ids = zeros(num_rois,1);

IdNames={};

for kk=1:1:num_rois
    IDs(kk) = str2num(xml1.labelset.label{kk}.Attributes.id);
    IDNames{kk} = xml1.labelset.label{kk}.Attributes.fullname;
end

LookUpTable.IDs=IDs;
LookUpTable.IDNames=IDNames;

save('BCI-DNI_brain_grayordinate_labels.mat', 'LookUpTable', 'labels');






