function tNLMPDFGOrdfMRI(GOrdInFile,GOrdOutFile,config)

load(GOrdInFile);
option = tNLMPdf();
option.FPR = str2double(config.fpr);
if ~strcmpi(config.memory,'auto')
    config.memory=str2double(config.memory);
end
option.memoryLimit = config.memory;
option.SCBFile = config.scbPath;
dtseries=tNLMPdf(dtseries, option);
save(GOrdOutFile,'dtseries');
