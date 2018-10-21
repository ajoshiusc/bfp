function tNLMPDFGOrdfMRI(GOrdInFile,GOrdOutFile,fpr,memory,MultiThreading)

load(GOrdInFile);
option = tNLMGPDF();
option.FPR = str2double(fpr);
if ~strcmpi(memory,'auto')
    memory=str2double(memory);
end

if ischar(MultiThreading)
    MultiThreading = str2double(MultiThreading);
end

if MultiThreading == 0
    option.numCPU=1;
else
    option.numCPU='auto';
end

option.memoryLimit = config.memory;
option.SCBFile = config.scbPath;
dtseries=tNLMGPDF(dtseries, option);
save(GOrdOutFile,'dtseries');
