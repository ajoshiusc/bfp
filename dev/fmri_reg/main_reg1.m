%||AUM||
%||Shree Ganeshaya Namaha||

clc;clear all;close all;
addpath(genpath('/home/ajoshi/coding_ground/bfp/dev'));
addpath(genpath('/home/ajoshi/coding_ground/bfp/src'));
%addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));

BFP_PATH='/home/ajoshi/coding_ground/bfp';
BrainSuitePath='/home/ajoshi/BrainSuite17a/svreg';

sub1='/deneb_disk/studyforrest_bfp/sub-03/func/sub-03_ses-movie_task-movie_run-3_bold.32k.GOrd.mat';
sub2='/deneb_disk/studyforrest_bfp/sub-02/func/sub-02_ses-movie_task-movie_run-3_bold.32k.GOrd.mat';

h=tic;
[wsub,origmap,newmap]=fmri_demons(sub1,sub2,BFP_PATH,'left');
toc(h)

figure;
patch('faces',sl.faces,'vertices',[xmap,ymap],'facevertexcdata',sqrt((xmap2'-xmap).^2+(ymap2'-ymap).^2),'edgecolor','k','facecolor','interp');
axis equal;axis off;camlight;material dull;

figure;
patch('faces',sl.faces,'vertices',[xmap2',ymap2'],'facevertexcdata',sqrt((xmap2'-xmap).^2+(ymap2'-ymap).^2),'edgecolor','k','facecolor','interp');
axis equal;axis off;camlight;material dull;

figure;plot(costiter);
%save temp1

%%
a=load('/home/ajoshi/coding_ground/bfp/supp_data/bci_grayordinates_surf_ind.mat');
atl=readdfs('/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.left.mid.cortex.dfs');
sl.labels=atl.labels(a.ind_left(ind));
writedfs('out1.dfs',sl);
recolor_by_label('out1.dfs','/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain');
sl=readdfs('out1.dfs');

slsm=smooth_cortex_fast(sl,.1,600);

figure;
patch('faces',sl.faces,'vertices',slsm.vertices,'facevertexcdata',sl.vcolor,'edgecolor','none','facecolor','interp');
axis equal;axis off;camlight;material dull;
slw=sl;
slw.labels=griddata(xmap,ymap,sl.labels,xmap2',ymap2','nearest');
writedfs('out1.dfs',slw);
recolor_by_label('out1.dfs','/home/ajoshi/coding_ground/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain');
slw=readdfs('out1.dfs');

figure;
patch('faces',slw.faces,'vertices',slsm.vertices,'facevertexcdata',slw.vcolor,'edgecolor','none','facecolor','interp');
axis equal;axis off;camlight;material dull;


