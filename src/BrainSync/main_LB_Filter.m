
% This is a sample main file that reads two grayordinate fmri files generated by bfp and outputs the synced files.

clc;clear all;close all;
addpath(genpath('/home/ajoshi/coding_ground/svreg'));
addpath(genpath('../'));

sl=readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kleft.dfs');
sr=readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kright.dfs');
msk=load('/big_disk/ajoshi/HCP100-fMRI-NLM/HCP100-fMRI-NLM/reference/100307.LR_mask.mat');
cSZ=length(msk.LR_flag);
ll=dir('/deneb_disk/Beijing_Zhang_bfp/*.mat');
load LBO_right
LBOr=LBO;
load LBO_left
ind=[1:cSZ];
ids={'sn8133','sn4055','tr4277','sn7915','sn5895','sn7602','sn6012','tr3170','sn6594','sn7256'};

for jj=1:length(ids)%ll)
%    g1=load(['/deneb_disk/Beijing_Zhang_bfp/',ll(jj).name]);
    subid=ids{jj};
    fname=['/deneb_disk/from_Todd_Constable_Epilepsy_Processed/',subid,'/func/',subid,'_rest_bold.32k.GOrd.mat'];
    if ~exist(['/deneb_disk/from_Todd_Constable_Epilepsy_Processed/',subid,'/func/',subid,'_rest_bold.32k.GOrd.mat'],'file')
        continue;
    end
    
    g1=load(fname);

    cSZ=length(msk.LR_flag);
    dtseries=g1.dtseries;
    g1.dtseries=g1.dtseries(ind,:);
    g1o=g1;    
    g1.dtseries=g1.dtseries(msk.LR_flag,:);
    
    dtseries(ind(msk.LR_flag),:) = laplaceBeltramiSmooth(LBO, g1.dtseries, 40);

    g1.dtseries=g1o.dtseries(msk.LR_flag==0,:);    
    dtseries(ind(msk.LR_flag==0),:) = laplaceBeltramiSmooth(LBOr, g1.dtseries, 40);
    save(['/deneb_disk/from_Todd_Constable_Epilepsy_Processed/',subid,'/func/',subid,'_rest_bold.32k.GOrd_LB40.mat'],'dtseries');
    %save(['/deneb_disk/Beijing_Zhang_bfp/',ll(jj).name(1:end-4),'_LB40.mat'],'dtseries');
    jj
end
