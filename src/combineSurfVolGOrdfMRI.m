function combineSurfVolGOrdfMRI(GOrdSurfIndFile,GOrdVolFile,GOrdFile)

load(GOrdSurfIndFile);
load(GOrdVolFile);
lsz=size(left_GO_fMRI,1);
rsz=size(right_GO_fMRI,1);
GO_fmri=Vol_GO_fmri;
GO_fmri(1:lsz,:)=left_GO_fMRI;
GO_fmri(1+lsz:lsz+rsz,:)=right_GO_fMRI;

v=make_nii(GO_fmri);

save_nii(v,GOrdFile);

