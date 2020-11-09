clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/ImagePTE1/ajoshi/code_farm/svreg/src'));
addpath(genpath('/ImagePTE1/ajoshi/code_farm/svreg/3rdParty'));

load('../supp_data/USCBrain_grayordinate_labels.mat');
roilist=[1147,225, 373,247,186,187,188,189];%, 345,351];


s{1}=readdfs('../supp_data/bci32kleft.dfs');
%s{2}=readdfs('../supp_data/bci32kright.dfs');

s1=combine_surf(s);
s1= smooth_cortex_fast(s1,.1,1000);

%writedfs('../supp_data/bci32k.dfs',s1);

labels=mod(labels,1000);
roilist=mod(roilist,1000);
ia = ismember(labels,roilist);
labels = double(labels).*double(ia);
s1.labels = labels(1:length(s1.vertices));

writedfs('tmpbci32k.dfs',s1);

recolor_by_label('tmpbci32k.dfs','/ImagePTE1/ajoshi/code_farm/svreg/USCBrain','/ImagePTE1/ajoshi/code_farm/svreg/USCBrain/brainsuite_labeldescription.xml');
