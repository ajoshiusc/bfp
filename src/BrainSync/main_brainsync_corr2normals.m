
% This is a sample main file that reads two grayordinate fmri files generated by bfp and outputs the synced files.

clc;clear all;close all;
addpath(genpath('/home/ajoshi/coding_ground/svreg'));
addpath(genpath('../'));
load LBO_right
LBOr=LBO;
load LBO_left
ids={'sn8133','sn4055','tr4277','sn7915','sn5895','sn7602','sn6012','tr3170','sn6594','sn7256','sub05267','sub06880'};
for jj=1:length(ids)
    subid=ids{jj};
    fname=['/deneb_disk/from_Todd_Constable_Epilepsy_Processed/',subid,'/func/',subid,'_rest_bold.32k.GOrd.mat'];
    if ~exist(fname,'file')
        fname=['/deneb_disk/Beijing_Zhang_bfp/',subid,'_rest_bold.32k.GOrd.mat'];
        if ~exist(fname,'file')
            continue;
        end
    end
    g2=load(fname);
    %g2=load('/deneb_disk/Beijing_Zhang_bfp/sub04050_rest_bold.32k.GOrd.mat');
    sl=readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kleft.dfs');
    sl=smooth_cortex_fast(sl,.1,1500);
    sr=readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kright.dfs');
    sr=smooth_cortex_fast(sr,.1,1500);
    msk=load('/big_disk/ajoshi/HCP100-fMRI-NLM/HCP100-fMRI-NLM/reference/100307.LR_mask.mat');
    
    
    ll=dir('/deneb_disk/Beijing_Zhang_bfp/*GOrd.mat');
    rho_after=0;g2oo=g2;
    rho_afterL=0;rho_afterR=0;
    for kk=1:length(ll)
        g1=load(['/deneb_disk/Beijing_Zhang_bfp/',ll(kk).name]);
        g1.dtseries=g1.dtseries(:,1:size(g2.dtseries,2));
        g1.dtseries=normalizeData(g1.dtseries')';
        g2=g2oo;
        
        cSZ=length(msk.LR_flag);
        g1.dtseries=g1.dtseries(1:cSZ,:);
        g2.dtseries=g2.dtseries(1:cSZ,:);
        
        g1o=g1;g2o=g2;
        
        g1.dtseries=g1.dtseries(msk.LR_flag,:);
        g2.dtseries=g2.dtseries(msk.LR_flag,:);
        
        
        g1.dtseries = laplaceBeltramiSmooth(LBO, g1.dtseries, 80);
        g2.dtseries = laplaceBeltramiSmooth(LBO, g2.dtseries, 80);
        
        g1.dtseries=normalizeData(g1.dtseries')';
        g2.dtseries=normalizeData(g2.dtseries')';
        rho_before=sum(g1.dtseries.*g2.dtseries,2);
        
        dtseries_sync=brainSync(g1.dtseries',g2.dtseries')';
        rho_after=sum(g1.dtseries.*dtseries_sync,2);
        rho_afterL=rho_afterL+rho_after;
        fprintf('Left Correlation before=%g, after=%g\n',mean(rho_before), mean(rho_after));
        
        %Right Hemisphere
        g1.dtseries=g1o.dtseries(msk.LR_flag==0,:);
        g2.dtseries=g2o.dtseries(msk.LR_flag==0,:);
        
        g1.dtseries = laplaceBeltramiSmooth(LBOr, g1.dtseries, 80);
        g2.dtseries = laplaceBeltramiSmooth(LBOr, g2.dtseries, 80);
        
        g1.dtseries=normalizeData(g1.dtseries')';
        g2.dtseries=normalizeData(g2.dtseries')';
        
        rho_before=sum(g1.dtseries.*g2.dtseries,2);
        
        dtseries_sync=brainSync(g1.dtseries',g2.dtseries')';
        rho_after=sum(g1.dtseries.*dtseries_sync,2);
        rho_afterR=rho_afterR+rho_after;
        fprintf('%d,%d Right Correlation before=%g, after=%g\n',jj,kk,mean(rho_before), mean(rho_after));
        
    end
    
    rho_afterL=rho_afterL/length(ll);
    rho_afterR=rho_afterR/length(ll);

    save(['corr_',subid,'_LB80.mat'],'rho_afterR','rho_afterL');

    h=figure;
    patch('faces',sr.faces,'vertices',sr.vertices,'facevertexcdata',rho_afterR,'facecolor','interp','edgecolor','none');
    view(90,0);axis equal;axis off;camlight;material dull;caxis([0,1]);colormap jet;
    saveas(h,['right_corr_',subid,'_LB80_1.png']);
    view(-90,0);axis equal;axis off;camlight;material dull;caxis([0,1]);colormap jet;
    saveas(h,['right_corr_',subid,'_LB80_2.png']);
    
    h=figure;
    patch('faces',sl.faces,'vertices',sl.vertices,'facevertexcdata',rho_afterL,'facecolor','interp','edgecolor','none');
    view(-90,0);axis equal;axis off;camlight;material dull;caxis([0,1]);colormap jet;
    saveas(h,['left_corr_',subid,'_LB80_1.png']);
    view(90,0);axis equal;axis off;camlight;material dull;caxis([0,1]);colormap jet;
    saveas(h,['left_corr_',subid,'_LB80_2.png']);
    
    figure;
    patch('faces',sl.faces,'vertices',sl.vertices,'facevertexcdata',rho_afterL-rho_afterR,'facecolor','interp','edgecolor','none');
    view(-90,0);axis equal;axis off;camlight;material dull;colormap jet;
    
end


