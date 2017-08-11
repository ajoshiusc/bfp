function tNLMPDFGOrdfMRI(GOrdInFile,GOrdOutFile)

load(GOrdInFile);
dtseries=tNLMPdf(dtseries);
save(v,GOrdOutFile);
