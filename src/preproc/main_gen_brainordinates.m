%% AUM %%
%Shree Ganeshaya Namaha

%function gen_brainordinates(subname,studydir)
clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/ImagePTE1/ajoshi/code_farm/bfp/dev/gifti-1.6'));
addpath(genpath('/ImagePTE1/ajoshi/code_farm/svreg/src'));
addpath(genpath('/ImagePTE1/ajoshi/code_farm/svreg/3rdParty'));

BSTDir = '/home/ajoshi/BrainSuite19b';
bfp_out ='/ImagePTE1/ajoshi/for_abhijit/bfp_out';
subid = 'sub-01';
sessionid='task-rest';

gen_brainordinates(BSTDir, bfp_out, subid, sessionid);

