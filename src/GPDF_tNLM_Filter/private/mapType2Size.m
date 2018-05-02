%
% map = mapType2Size(inByte)
% 
% Description:
%     generate a map from type of variable to Byte
% 
% Input:
%     inByte - flag indicating whether in Byte or Bit
% 
% Output:
%     map - the map
% 
% Copyright:
%     2017-2018 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     1.0.2
% Date:
%     2017/08/16
%

function map = mapType2Size(inByte)

    if ~exist('inByte', 'var') || isempty(inByte)
        inByte = true;
    end

    allTypes = {'int8', 'int16', 'int32', 'int64', ...
                'uint8', 'uint16', 'uint32' 'uint64', ...
                'single', 'double', ...
                'char', ...
                'logical'};
    
    allSize = {8, 16, 32, 64, ...
               8, 16, 32, 64, ...
               32, 64, ...
               16, ...
               8};
    
    if inByte
        allSize = cellfun(@(x) x/8, allSize, 'UniformOutput', false);
    end
    
    map = containers.Map(allTypes, allSize);
    
end