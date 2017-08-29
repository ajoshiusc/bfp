clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/bfp/dev/gifti-1.6'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/3rdParty'));
addpath(genpath('/big_disk/ajoshi/freesurfer/matlab'));

%% Resample from one sphere to the other (Right)
g_32k = gifti('fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii');
[bci.vertices, bci.faces] = freesurfer_read_surf('/big_disk/ajoshi/data/BCI_DNI_Atlas/surf/rh.sphere.reg');
tic
[ind_right,d] = dsearchn(bci.vertices,g_32k.vertices);
toc
g.vertices=bci.vertices(ind_right,:);
g.faces=g_32k.faces;

%% Resample to 32k Surfaces (Right)
[bci_inner_fs.vertices, bci_inner_fs.faces] = freesurfer_read_surf('/big_disk/ajoshi/data/BCI_DNI_Atlas/surf/rh.orig');
bci_BST=readdfs('/home/ajoshi/BrainSuite17a/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.right.inner.cortex.dfs');
bci_BST.vertices(:, 1) = bci_BST.vertices(:, 1)-96*0.8;
bci_BST.vertices(:, 2) = bci_BST.vertices(:, 2)-192*0.546875;
bci_BST.vertices(:, 3) = bci_BST.vertices(:, 3)-192*0.546875;

tic
[ind_right,d] = dsearchn(bci_BST.vertices,bci_inner_fs.vertices(ind_right,:));
toc
bci_BST=readdfs('/home/ajoshi/BrainSuite17a/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.right.mid.cortex.dfs');

bci32kright.faces=double(g.faces);
bci32kright.vertices=bci_BST.vertices(ind_right,:);
bci32kright=smooth_cortex_fast(bci32kright,.1,75);

view_patch(bci32kright);
writedfs('../supp_data/bci32kright.dfs',bci32kright);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resample from one sphere to the other (Left)
g_32k = gifti('fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii');
[bci.vertices, bci.faces] = freesurfer_read_surf('/big_disk/ajoshi/data/BCI_DNI_Atlas/surf/lh.sphere.reg');
tic
[ind_left,d] = dsearchn(bci.vertices,g_32k.vertices);
toc
g.vertices=bci.vertices(ind_left,:);
g.faces=g_32k.faces;

%% Resample to 32k Surfaces (Left)
[bci_inner_fs.vertices, bci_inner_fs.faces] = freesurfer_read_surf('/big_disk/ajoshi/data/BCI_DNI_Atlas/surf/lh.orig');
bci_BST=readdfs('/home/ajoshi/BrainSuite17a/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.left.inner.cortex.dfs');
bci_BST.vertices(:, 1) = bci_BST.vertices(:, 1)-96*0.8;
bci_BST.vertices(:, 2) = bci_BST.vertices(:, 2)-192*0.546875;
bci_BST.vertices(:, 3) = bci_BST.vertices(:, 3)-192*0.546875;

tic
[ind_left,d] = dsearchn(bci_BST.vertices,bci_inner_fs.vertices(ind_left,:));
toc
bci_BST=readdfs('/home/ajoshi/BrainSuite17a/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.left.mid.cortex.dfs');

bci32kleft.faces=double(g.faces);
bci32kleft.vertices=bci_BST.vertices(ind_left,:);

bci32kleft=smooth_cortex_fast(bci32kleft,.1,75);
view_patch(bci32kleft);

writedfs('../supp_data/bci32kleft.dfs',bci32kleft);

save('../supp_data/bci_grayordinates_surf_ind.mat','ind_left','ind_right');
