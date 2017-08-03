function generateSurfGOrdfMRI(GOrdIndFile,subbasename)

load(GOrdIndFile);

load([subbasename,'_fmri2surf.mat']);

left_GO_fMRI=datal_atlas(ind_left,:);
right_GO_fMRI=datar_atlas(ind_right,:);

save([subbasename,'_fmri2surf.mat'],'left_GO_fMRI','right_GO_fMRI');

