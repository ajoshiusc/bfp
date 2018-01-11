%||AUM||
%||Shree Ganeshaya Namaha||

clc;clear all;close all;
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/3rdParty'));

addpath(genpath('/home/ajoshi/coding_ground/bfp/dev'));
addpath(genpath('/home/ajoshi/coding_ground/bfp/src'));

BFPPATH='/home/ajoshi/coding_ground/bfp';
BrainSuitePath='/home/ajoshi/BrainSuite17a/svreg';

sub1='/deneb_disk/studyforrest_bfp/sub-03/func/sub-03_ses-movie_task-movie_run-3_bold.32k.GOrd.mat';
sub2='/deneb_disk/studyforrest_bfp/sub-02/func/sub-02_ses-movie_task-movie_run-3_bold.32k.GOrd.mat';

h=tic;
[wsub,origmap,newmap,s,costiter,ind]=fmri_demons(sub1,sub2,BFPPATH,'left');
toc(h)

save('aaj1000.mat','wsub','origmap','newmap','s','costiter','ind');

[~,C1]=vertices_connectivity_fast(s);
[s,A1]=smooth_cortex_fast(s,.8,50);
s.attributes=curvature_cortex_fast(s,50,0,C1);

figure;
patch('faces',s.faces,'vertices',origmap,'facevertexcdata',sqrt(sum((origmap-newmap).^2,2)),'edgecolor','none','facecolor','interp');
axis equal;axis off;camlight;material dull;

figure;
patch('faces',s.faces,'vertices',newmap,'facevertexcdata',sqrt(sum((origmap-newmap).^2,2)),'edgecolor','k','facecolor','interp');
axis equal;axis off;camlight;material dull;

figure;plot(costiter);
%save temp1

%%
a=load('/home/ajoshi/coding_ground/bfp/supp_data/bci_grayordinates_surf_ind.mat');
atl=readdfs('/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.left.mid.cortex.dfs');
s.labels=atl.labels(a.ind_left(ind));
writedfs('out1.dfs',s);
recolor_by_label('out1.dfs','/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain');
s=readdfs('out1.dfs');

sm=smooth_cortex_fast(s,.1,600);

figure;
patch('faces',s.faces,'vertices',sm.vertices,'facevertexcdata',s.vcolor,'edgecolor','none','facecolor','interp');
axis equal;axis off;view(-90,0);camlight;material dull;
sw=s;
sw.labels=griddata(origmap(:,1),origmap(:,2),s.labels,newmap(:,1),newmap(:,2),'nearest');
writedfs('out1.dfs',sw);
recolor_by_label('out1.dfs','/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain');
sw=readdfs('out1.dfs');

figure;
patch('faces',sw.faces,'vertices',sm.vertices,'facevertexcdata',sw.vcolor,'edgecolor','none','facecolor','interp');
axis equal;axis off;view(-90,0);camlight;material dull;


