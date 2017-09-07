function generateSurfGOrdfMRI(GOrdIndFile,fmri2surfFile,GOrdsurfFile)

load(GOrdIndFile);
load(fmri2surfFile);

left_GO_fMRI=datal_atlas(ind_left,:);
right_GO_fMRI=datar_atlas(ind_right,:);

save(GOrdsurfFile,'left_GO_fMRI','right_GO_fMRI');

