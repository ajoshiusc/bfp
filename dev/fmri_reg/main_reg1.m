%||AUM||
%||Shree Ganeshaya Namaha||

clc;clear all;close all;
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/3rdParty'));

addpath(genpath('/home/ajoshi/coding_ground/bfp/dev'));
addpath(genpath('/home/ajoshi/coding_ground/bfp/src'));

BFPPATH='/home/ajoshi/coding_ground/bfp';
BrainSuitePath='/home/ajoshi/BrainSuite17a/svreg';

for sc=1:2
    sub1=sprintf('/big_disk/ajoshi/HCP5/100307/MNINonLinear/Results/rfMRI_REST%d_LR/rfMRI_REST%d_LR_Atlas_hp2000_clean.dtseries.nii',sc,sc);
    sub2=sprintf('/big_disk/ajoshi/HCP5/103414/MNINonLinear/Results/rfMRI_REST%d_LR/rfMRI_REST%d_LR_Atlas_hp2000_clean.dtseries.nii',sc,sc);
    sub1=ft_read_cifti(sub1);
    sub2=ft_read_cifti(sub2);
    
    % a=load('/big_disk/ajoshi/HCP5-fMRI-NLM/100307/100307.rfMRI_REST1_LR.reduce3.ftdata.NLM_11N_hvar_25.mat');
    % sub1.dtseries=a.ftdata_NLM;
    % a=load('/big_disk/ajoshi/HCP5-fMRI-NLM/110411/110411.rfMRI_REST1_LR.reduce3.ftdata.NLM_11N_hvar_25.mat');
    % sub2.dtseries=a.ftdata_NLM;
    
    %sub1.dtseries=sub1.dtseries(:,[1:100]);
    %sub2.dtseries=sub2.dtseries(:,[1:100]);
    h=tic;
    [wsub,origmap,newmap,s,costiter,ind]=fmri_demons(sub1,sub2,BFPPATH,'left');
    t1=toc(h)
    
    save(sprintf('aaj1200_%d.mat',sc),'wsub','origmap','newmap','s','costiter','ind');
    
    [~,C1]=vertices_connectivity_fast(s);
    [s,A1]=smooth_cortex_fast(s,.8,50);
    s.attributes=curvature_cortex_fast(s,50,0,C1);
    
    h=figure;
    patch('faces',s.faces,'vertices',origmap,'facevertexcdata',sqrt(sum((origmap-newmap).^2,2)),'edgecolor','none','facecolor','interp');
    axis equal;axis off;material dull;axis tight;
    saveas(h,'sub_sqr_mesh.png');
    h=figure;
    patch('faces',s.faces,'vertices',newmap,'facevertexcdata',sqrt(sum((origmap-newmap).^2,2)),'edgecolor','y','facecolor','interp');
    axis equal;axis off;material dull;axis tight;
    saveas(h,'warped_sqr_mesh.png');
    
    figure;plot(costiter);
    %save temp1
    
    %%
    a=load('/home/ajoshi/coding_ground/bfp/supp_data/bci_grayordinates_surf_ind.mat');
    atl=readdfs('/home/ajoshi/coding_ground/svreg/USCBrain/BCI-DNI_brain.left.mid.cortex.dfs');
    s.labels=atl.labels(a.ind_left(ind));
    writedfs('out1.dfs',s);
    recolor_by_label('out1.dfs','/home/ajoshi/coding_ground/svreg/USCBrain/BCI-DNI_brain');
    s=readdfs('out1.dfs');
    
    sm=smooth_cortex_fast(s,.1,600);
    sm.faces=sm.faces(:,[2,1,3]);
    sf=figure;
    patch('faces',sm.faces,'vertices',sm.vertices,'facevertexcdata',s.vcolor,'edgecolor','none','facecolor','interp','backfacelighting','unlit');
    axis equal;axis off;view(-90,0);camlight;material dull;axis tight;
    saveas(sf,sprintf('sub_left1_%d.png',sc));
    view(90,0);camlight;
    saveas(sf,sprintf('sub_left2_%d.png',sc));
    
    sw=s;
    sw.labels=griddata(origmap(:,1),origmap(:,2),s.labels,newmap(:,1),newmap(:,2),'nearest');
    writedfs('out1.dfs',sw);
    recolor_by_label('out1.dfs','/home/ajoshi/coding_ground/svreg/USCBrain/BCI-DNI_brain');
    sw=readdfs('out1.dfs');
    
    w=figure;
    patch('faces',sm.faces,'vertices',sm.vertices,'facevertexcdata',sw.vcolor,'edgecolor','none','facecolor','interp','backfacelighting','unlit');
    axis equal;axis off;view(-90,0);camlight;material dull;axis tight;
    saveas(w,sprintf('warped_left1_%d.png',sc));
    view(90,0);camlight;
    saveas(w,sprintf('warped_left2_%d.png',sc));
    
    %%
    % This is not meaningful since the deformation is on the square
    w=figure;
    patch('faces',sm.faces,'vertices',sm.vertices,'facevertexcdata',sqrt(sum((origmap-newmap).^2,2)),'edgecolor','none','facecolor','interp','backfacelighting','unlit');
    axis equal;axis off;view(-90,0);camlight;material dull;axis tight;
    saveas(w,'def_left1.png');
    view(90,0);camlight;
    saveas(w,'def_left2.png');
    
    
    sq=figure;
    patch('faces',s.faces,'vertices',origmap,'facevertexcdata',s.vcolor,'edgecolor','none','facecolor','interp');
    axis equal;axis off;camlight;material dull;axis tight;
    saveas(sq,'flat.png');
    sq=figure;
    patch('faces',s.faces,'vertices',newmap,'facevertexcdata',sw.vcolor,'edgecolor','none','facecolor','interp');
    axis equal;axis off;camlight;material dull;
    saveas(sq,'warped_flat.png');
    
end
