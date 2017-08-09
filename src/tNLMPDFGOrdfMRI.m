function tNLMPDFGOrdfMRI(GOrdInFile,GOrdOutFile)

v=load_nii(v,GOrdInFile);
v.img=tNLMPDF(v.img);
save_nii(v,GOrdOutFile);
