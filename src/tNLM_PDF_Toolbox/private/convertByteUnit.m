%
% [dataOut, unitOut] = convertByteUnit(dataIn, from, to)
% 
% Description:
%     convert unit of Byte from one to another
% 
% Input:
%     dataIn - the number of Bytes in "from" unit
%     from - unit for current number
%     to - unit to convert to
% 
% Output:
%     dataOut - converted data
%     unitOut - unit converted to, useful if "to" is set to auto
% 
% Copyright:
%     2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     1.0.2
% Date:
%     2017/08/16
%

function [dataOut, unitOut] = convertByteUnit(dataIn, from, to)
    
    units = getByteUnits();
    
    if strcmp(to, 'auto')
        dataInByte = convertByteUnit(dataIn, from, 'B');
        sc = floor(log(dataInByte)/log(1024));
        to = units{sc+1};
    end
    
    if ~(ismember(from, units) && ismember(to, units))
        disp(units);
        error('unit has to be one of above');
    end
    
    idxF = find(ismember(units, from));
    idxT = find(ismember(units, to));
    idxD = idxF - idxT;
    
    dataOut = dataIn .* (1024^idxD);
    unitOut = to;
    
end