
% The purpose of this file is to generate a file that gives a vector of
% grayordinate labels, its id lookup tables
clc;clear all;close all;
addpath(genpath('/ImagePTE1/ajoshi/code_farm/bfp/dev/gifti-1.6'));
addpath(genpath('/ImagePTE1/ajoshi/code_farm/svreg/src'));
addpath(genpath('/ImagePTE1/ajoshi/code_farm/svreg/3rdParty'));

% Read surface indices
load('../supp_data/bci_grayordinates_surf_ind.mat')
load('../supp_data/HCP_32k_Label.mat')
hcp_isnan=isnan(brainstructure(1:2*32492));

l= readdfs('/ImagePTE1/ajoshi/code_farm/svreg/USCBrain/USCBrain.left.mid.cortex.dfs');

r= readdfs('/ImagePTE1/ajoshi/code_farm/svreg/USCBrain/USCBrain.right.mid.cortex.dfs');

surf_labs=[l.labels(ind_left);r.labels(ind_right)];

surf_labs(hcp_isnan)=0;

% Read indices for volume and read labels from the corresponding volume

load('../supp_data/bci_grayordinates_vol_ind.mat')
bci_vol_ind(isnan(bci_vol_ind))=1;
vl = load_nii('/ImagePTE1/ajoshi/code_farm/svreg/USCBrain/USCBrain.label.nii.gz');

vol_labs=vl.img(bci_vol_ind);

% Concatenate surface and volume labels
labels=[surf_labs;vol_labs(1+length(surf_labs):end)];

% Read XML and create lookup table
xml1 = xml2struct('/ImagePTE1/ajoshi/code_farm/svreg/USCBrain/brainsuite_labeldescription.xml');

num_rois = length(xml1.labelset.label);
Ids = zeros(num_rois,1);

IdNames={};

for kk=1:1:num_rois
    IDs(kk) = str2num(xml1.labelset.label{kk}.Attributes.id);
    IDNames{kk} = xml1.labelset.label{kk}.Attributes.fullname;
end

LookUpTable.IDs=IDs;
LookUpTable.IDNames=IDNames;

save('USCBrain_grayordinate_labels.mat', 'LookUpTable', 'labels');







