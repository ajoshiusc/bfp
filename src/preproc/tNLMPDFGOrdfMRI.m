function tNLMPDFGOrdfMRI(GOrdInFile,GOrdOutFile)

v=load_nii(GOrdInFile);
v.img=tNLMPdf(v.img);
save_nii(v,GOrdOutFile);
