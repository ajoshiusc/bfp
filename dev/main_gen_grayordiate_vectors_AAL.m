
% The purpose of this file is to generate a file that gives a vector of
% grayordinate labels, its id lookup tables
clc;clear all;close all;restoredefaultpath;
addpath(genpath('/ImagePTE1/ajoshi/code_farm/bfp/dev/gifti-1.6'));
addpath(genpath('/ImagePTE1/ajoshi/code_farm/bfp/src'));
%addpath(genpath('/ImagePTE1/ajoshi/code_farm/svreg/3rdParty'));

% Read surface indices
load('../supp_data/colin27_grayordinates_surf_ind.mat')
load('../supp_data/HCP_32k_Label.mat')
hcp_isnan=isnan(brainstructure(1:2*32492));

l= readdfs('/ImagePTE1/ajoshi/code_farm/hybridatlas/src/colin27-AAL.left.mid.cortex.dfs');

r= readdfs('/ImagePTE1/ajoshi/code_farm/hybridatlas/src/colin27-AAL.right.mid.cortex.dfs');

surf_labs=[l.labels(ind_left);r.labels(ind_right)];

surf_labs(hcp_isnan)=0;

% Read indices for volume and read labels from the corresponding volume

load('../supp_data/colin27_grayordinates_vol_ind.mat')
col_vol_ind(isnan(col_vol_ind))=1;
vl = load_nii('/ImagePTE1/ajoshi/code_farm/hybridatlas/aal_for_SPM8/aal_for_SPM8/ROI_MNI_V4_svreg_1mm.nii.gz');

vol_labs=vl.img(col_vol_ind);

% Concatenate surface and volume labels
labels=[surf_labs;vol_labs(1+length(surf_labs):end)];

% Read XML and create lookup table
xml1 = readtable('/ImagePTE1/ajoshi/code_farm/hybridatlas/aal_for_SPM8/aal_for_SPM8/ROI_MNI_V4.txt');

num_rois = length(xml1{:,1});
Ids = zeros(num_rois,1);

IdNames={};

for kk=1:1:num_rois
    IDs(kk) = xml1{kk,3};
    IDNames{kk} = xml1{kk,2}{:};
end

LookUpTable.IDs=IDs;
LookUpTable.IDNames=IDNames;

save('AAL_grayordinate_labels.mat', 'LookUpTable', 'labels');

