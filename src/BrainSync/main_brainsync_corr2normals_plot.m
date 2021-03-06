
% This is a sample main file that reads two grayordinate fmri files generated by bfp and outputs the synced files.

clc;clear all;close all;
addpath(genpath('/home/ajoshi/coding_ground/svreg'));
addpath(genpath('../'));
load LBO_right
LBOr=LBO;
load LBO_left
ids={'sub18604','sn4055','sub08224','sub05676','sn8133','tr4277','sn7915','sn5895','sn7602','sn6012','tr3170','sn6594','sn7256','sub05267','sub06880'};

    sl=readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kleft.dfs');
    sl=smooth_cortex_fast(sl,.1,500);
    sr=readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kright.dfs');
    sr=smooth_cortex_fast(sr,.1,500);
hcpref=readdfs('/big_disk/ajoshi/HCP_data/reference/196750.aparc.a2009s.32k_fs.left.dfs.gz');
figure;
patch('faces',sl.faces,'vertices',sl.vertices,'facevertexcdata',1.0*(hcpref.labels==0),'edgecolor','none','facecolor','interp');
axis equal;axis off;
CCind=find(hcpref.labels==0);

for jj=1:length(ids)
    subid=ids{jj};
    fname=['/deneb_disk/from_Todd_Constable_Epilepsy_Processed/',subid,'/func/',subid,'_rest_bold.32k.GOrd.mat'];
    if ~exist(fname,'file')
        fname=['/deneb_disk/Beijing_Zhang_bfp/',subid,'_rest_bold.32k.GOrd.mat'];
        if ~exist(fname,'file')
            continue;
        end
    end

    load(['corr_',subid,'_LB80.mat'],'rho_afterR','rho_afterL');
    rho_afterR(CCind)=0;
    h=figure;
    patch('faces',sr.faces,'vertices',sr.vertices,'facevertexcdata',rho_afterR,'facecolor','interp','edgecolor','none');
    view(90,0);axis equal;axis off;camlight;material dull;caxis([0.5,1]);colormap jet;
    saveas(h,['right_corr_',subid,'_LB80_1_new.png']);
    view(-90,0);axis equal;axis off;camlight;material dull;caxis([0.5,1]);colormap jet;
    saveas(h,['right_corr_',subid,'_LB80_2_new.png']);

    rho_afterR(CCind)=1;
    h=figure;
    patch('faces',sr.faces,'vertices',sr.vertices,'facevertexcdata',1.0*(rho_afterR>0.8),'facecolor','interp','edgecolor','none');
    view(90,0);axis equal;axis off;camlight;material dull;caxis([0.5,1]);colormap jet;
    saveas(h,['right_corr_',subid,'_LB80_1_new_pval.png']);
    view(-90,0);axis equal;axis off;camlight;material dull;caxis([0.5,1]);colormap jet;
    saveas(h,['right_corr_',subid,'_LB80_2_new_pval.png']);
    
    rho_afterL(CCind)=0;
    h=figure;
    patch('faces',sl.faces,'vertices',sl.vertices,'facevertexcdata',rho_afterL,'facecolor','interp','edgecolor','none');
    view(-90,0);axis equal;axis off;camlight;material dull;caxis([0.5,1]);colormap jet;
    saveas(h,['left_corr_',subid,'_LB80_1_new.png']);
    view(90,0);axis equal;axis off;camlight;material dull;caxis([0.5,1]);colormap jet;
    saveas(h,['left_corr_',subid,'_LB80_2_new.png']);

    rho_afterL(CCind)=1;
    h=figure;
    patch('faces',sl.faces,'vertices',sl.vertices,'facevertexcdata',1.0*(rho_afterL>0.8),'facecolor','interp','edgecolor','none');
    view(-90,0);axis equal;axis off;camlight;material dull;caxis([0.5,1]);colormap jet;
    saveas(h,['left_corr_',subid,'_LB80_1_new_pval.png']);
    view(90,0);axis equal;axis off;camlight;material dull;caxis([0.5,1]);colormap jet;
    saveas(h,['left_corr_',subid,'_LB80_2_new_pval.png']);
    
close all;    
end


