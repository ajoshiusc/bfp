function combineSurfVolGOrdfMRI(GOrdSurfFile,GOrdVolFile,GOrdFile)

load(GOrdSurfFile);
load(GOrdVolFile);
lsz=size(left_GO_fMRI,1);
rsz=size(right_GO_fMRI,1);
GO_fmri=Vol_GO_fmri;
GO_fmri(1:lsz,:)=left_GO_fMRI;
GO_fmri(1+lsz:lsz+rsz,:)=right_GO_fMRI;

dtseries=GO_fmri;
save(GOrdFile,'dtseries');

