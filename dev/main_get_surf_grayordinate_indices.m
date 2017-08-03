clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/bfp/dev/gifti-1.6'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/3rdParty'));
addpath(genpath('/big_disk/ajoshi/freesurfer/matlab'));

bci_BST = readdfs('/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.right.inner.cortex.dfs');
bci_BST.vertices(:, 1) = bci_BST.vertices(:, 1)-96*0.8;
bci_BST.vertices(:, 2) = bci_BST.vertices(:, 2)-192*0.546875;
bci_BST.vertices(:, 3) = bci_BST.vertices(:, 3)-192*0.546875;

[bci.vertices, bci.faces] = freesurfer_read_surf('/big_disk/ajoshi/data/BCI_DNI_Atlas/surf/rh.white')
view_patch(bci)

%''' FS_BCI to FS BCI Sphere'''
[bci.vertices, bci.faces] = freesurfer_read_surf('/big_disk/ajoshi/data/BCI_DNI_Atlas/surf/rh.sphere.reg')

%    ''' FS BCI Sphere to ref FS Sphere'''
g_surf = gifti('/big_disk/ajoshi/HCP_data/reference.old/100307/MNINonLinear/Native/100307.R.sphere.reg.native.surf.gii')
s.vertices = g_surf.vertices
s.faces = g_surf.faces

g_surf = gifti('/big_disk/ajoshi/HCP_data/reference.old/100307/MNINonLinear/Native/100307.R.very_inflated.native.surf.gii')

gr = readdfs('/big_disk/ajoshi/HCP_data/reference.old/100307.aparc.a2009s.32k_fs.very_smooth.right.dfs.gz')
