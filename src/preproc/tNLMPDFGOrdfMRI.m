function tNLMPDFGOrdfMRI(GOrdInFile,GOrdOutFile,config)

load(GOrdInFile);
option = tNLMGPDF();
option.FPR = str2double(config.fpr);
if ~strcmpi(config.memory,'auto')
    config.memory=str2double(config.memory);
end

if config.MultiThreading == 0
    option.numCPU=1;
else
    option.numCPU='auto';
end

option.memoryLimit = config.memory;
option.SCBFile = config.scbPath;
dtseries=tNLMGPDF(dtseries, option);
save(GOrdOutFile,'dtseries');
