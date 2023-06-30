function [pathstr, baseName] = fileBaseName( fileName )
% Returns path string & the base name of the fileName (removing multiple
% extensions - see remove_extension.m)
%
% Usage:
%     [pathstr, baseName ] = file_base_name(fileName)
%     baseName = file_base_name(fileName)
%


[pathstr, baseName, ext] = fileparts(fileName);
baseName = remove_extension([baseName ext]);

if (nargout <= 1)
   pathstr = baseName;
end


end

