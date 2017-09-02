function tNLMPDFGOrdfMRI(GOrdInFile,GOrdOutFile)

load(GOrdInFile);
dtseries=tNLMPdf(dtseries);
save(GOrdOutFile,'dtseries');
