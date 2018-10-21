 %#function gcp

function tNLMPDFGOrdfMRI(GOrdInFile,GOrdOutFile,fpr,memory,MultiThreading, scbPath)

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

option.memoryLimit = memory;
option.SCBFile = scbPath;
dtseries=tNLMGPDF(dtseries, option);
save(GOrdOutFile,'dtseries');
