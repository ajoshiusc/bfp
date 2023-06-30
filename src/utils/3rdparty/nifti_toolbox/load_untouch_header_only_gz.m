function [header, ext, filetype, machine] = load_untouch_header_only_gz( file )
% wrapper around load_untouch_header_only to work with .gz files.
%
%   1. Extracts GNU zip files, if file is zipped.
%   2. Loads the data & returns it.
%   3. Deletes the unzipped file, if file is zipped.
%
%   file - is a string with relative / absolute path to .gz or .nii file
%
%   Requires NIFTI toolbox.

loc = strfind(file, '.gz');
len = length(file);

if (~isempty(loc)) && (len-loc(end) == 2)  % .gz extension found
   workdir = tempname;
   niiData = gunzip(file, workdir);
   [header, ext, filetype, machine] = load_untouch_header_only(char(niiData));
   delete(char(niiData));
   rmdir(workdir, 's')
else
  [header, ext, filetype, machine] = load_untouch_header_only(file);
end

end

