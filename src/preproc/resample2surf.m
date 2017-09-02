% This file resamples fmri to surface
% subbasename is for the structural T1 image
%fmri is the file name of the fmri
% 
function resample2surf(subbasename,fmri,fmri2surfFile)

subpath=fileparts(subbasename);

v=load_nii_BIG_Lab(fmri);
sl=readdfs([subbasename,'.left.mid.cortex.svreg.dfs']);
sr=readdfs([subbasename,'.right.mid.cortex.svreg.dfs']);
al=readdfs([subpath,'/atlas.left.mid.cortex.svreg.dfs']);
ar=readdfs([subpath,'/atlas.right.mid.cortex.svreg.dfs']);

datal=zeros(length(sl.vertices),size(v.img,4));
datar=zeros(length(sr.vertices),size(v.img,4));
datal_atlas=zeros(length(al.vertices),size(v.img,4));
datar_atlas=zeros(length(ar.vertices),size(v.img,4));

res=v.hdr.dime.pixdim(2:4);
vimg=double(v.img);
parfor j=1:size(v.img,4)
    datal(:,j)=interp3(vimg(:,:,:,j), sl.vertices(:,2)/res(2) + 1,sl.vertices(:,1)/res(1) + 1, sl.vertices(:,3)/res(3) + 1);
    datar(:,j)=interp3(vimg(:,:,:,j), sr.vertices(:,2)/res(2) + 1,sr.vertices(:,1)/res(1) + 1, sr.vertices(:,3)/res(3) + 1);
end

v1=var(datal,[],2);
v2=var(datar,[],2);
sl2=sl; v1=(v1==0); ind=(sum(v1(sl2.faces),2)>0); sl2.faces(ind,:)=[];[sl2,locl]=myclean_patch_cc(sl2);
sr2=sr; v2=(v2==0); ind=(sum(v2(sr2.faces),2)>0); sr2.faces(ind,:)=[];[sr2,locr]=myclean_patch_cc(sr2);

parfor j=1:size(v.img,4)
    datal_atlas(:,j)=map_data_flatmap(sl2,datal(locl,j),al);
    datar_atlas(:,j)=map_data_flatmap(sr2,datar(locr,j),ar);        
    fprintf('%d,',j);    
end

save(fmri2surfFile,'datal_atlas','datar_atlas','datal','datar');


