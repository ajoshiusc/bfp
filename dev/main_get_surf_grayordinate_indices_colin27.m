clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/ImagePTE1/ajoshi/code_farm/bfp/dev/gifti-1.6'));
addpath(genpath('/ImagePTE1/ajoshi/code_farm/svreg/src'));
addpath(genpath('/ImagePTE1/ajoshi/code_farm/svreg/3rdParty'));
addpath(genpath('/big_disk/ajoshi/freesurfer/matlab'));

%% Resample from one sphere to the other (Right)
g_32k = gifti('/big_disk/ajoshi/data/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii');
[col.vertices, col.faces] = freesurfer_read_surf('/ImagePTE1/ajoshi/mni_colin27_1998/mni_colin27_fs/surf/rh.sphere.reg');
tic
[ind_right,d] = dsearchn(col.vertices,g_32k.vertices);
toc
g.vertices=col.vertices(ind_right,:);
g.faces=g_32k.faces;

%% Resample to 32k Surfaces (Right)
[col_inner_fs.vertices, col_inner_fs.faces] = freesurfer_read_surf('/ImagePTE1/ajoshi/mni_colin27_1998/mni_colin27_fs/surf/rh.orig');
col_BST=readdfs('/ImagePTE1/ajoshi/mni_colin27_1998/mni_colin27_BST/colin27_t1_tal_lin.right.inner.cortex.svreg.dfs');
col_BST.vertices(:, 1) = col_BST.vertices(:, 1)-181.0/2.0;
col_BST.vertices(:, 2) = col_BST.vertices(:, 2)-217.0/2.0;
col_BST.vertices(:, 3) = col_BST.vertices(:, 3)-181.0/2.0;

tic
[ind_right,d] = dsearchn(col_BST.vertices,col_inner_fs.vertices(ind_right,:));
toc
col_BST=readdfs('/ImagePTE1/ajoshi/mni_colin27_1998/mni_colin27_BST/colin27_t1_tal_lin.right.mid.cortex.svreg.dfs');

col32kright.faces=double(g.faces);
col32kright.vertices=col_BST.vertices(ind_right,:);

writedfs('../supp_data/colin27_32kright_orig.dfs',col32kright);

col32kright=smooth_cortex_fast(col32kright,.1,75);

view_patch(col32kright);
writedfs('../supp_data/colin27_32kright.dfs',col32kright);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resample from one sphere to the other (Left)
g_32k = gifti('/big_disk/ajoshi/data/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii');
[col.vertices, col.faces] = freesurfer_read_surf('/ImagePTE1/ajoshi/mni_colin27_1998/mni_colin27_fs/surf/lh.sphere.reg');
tic
[ind_left,d] = dsearchn(col.vertices,g_32k.vertices);
toc
g.vertices=col.vertices(ind_left,:);
g.faces=g_32k.faces;

%% Resample to 32k Surfaces (Left)
[col_inner_fs.vertices, col_inner_fs.faces] = freesurfer_read_surf('/ImagePTE1/ajoshi/mni_colin27_1998/mni_colin27_fs/surf/lh.orig');
col_BST=readdfs('/ImagePTE1/ajoshi/mni_colin27_1998/mni_colin27_BST/colin27_t1_tal_lin.left.inner.cortex.svreg.dfs');
col_BST.vertices(:, 1) = col_BST.vertices(:, 1)-181.0/2.0;
col_BST.vertices(:, 2) = col_BST.vertices(:, 2)-217.0/2.0;
col_BST.vertices(:, 3) = col_BST.vertices(:, 3)-181.0/2.0;

tic
[ind_left,d] = dsearchn(col_BST.vertices,col_inner_fs.vertices(ind_left,:));
toc
col_BST=readdfs('/ImagePTE1/ajoshi/mni_colin27_1998/mni_colin27_BST/colin27_t1_tal_lin.left.mid.cortex.svreg.dfs');

col32kleft.faces=double(g.faces);
col32kleft.vertices=col_BST.vertices(ind_left,:);

writedfs('../supp_data/colin27_32kleft_orig.dfs',col32kleft);

col32kleft=smooth_cortex_fast(col32kleft,.1,75);
view_patch(col32kleft);

writedfs('../supp_data/colin27_32kleft.dfs',col32kleft);

save('../supp_data/colin27_grayordinates_surf_ind.mat','ind_left','ind_right');
