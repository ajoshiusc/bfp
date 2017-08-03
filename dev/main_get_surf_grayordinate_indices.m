clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/bfp/dev/gifti-1.6'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/3rdParty'));
addpath(genpath('/big_disk/ajoshi/freesurfer/matlab'));


g_32k = gifti('fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii');
[bci.vertices, bci.faces] = freesurfer_read_surf('/big_disk/ajoshi/data/BCI_DNI_Atlas/surf/rh.sphere.reg');
tic
[ind_right,d] = dsearchn(bci.vertices,g_32k.vertices);
toc

g_32k = gifti('fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii');
[bci.vertices, bci.faces] = freesurfer_read_surf('/big_disk/ajoshi/data/BCI_DNI_Atlas/surf/lh.sphere.reg');
tic
[ind_left,d] = dsearchn(bci.vertices,g_32k.vertices);
toc
% g_32k.vertices=bci.vertices(ind_left,:);
% view_patch(g_32k)
save('bci_grayordinates_surf_ind.mat','ind_left','ind_right');
