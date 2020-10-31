function [ data ] = load_untouch_eig_gz(file, workdir)
% Loads .eig file to nii structure. Each voxel has 12 elements which follows following eig
% structure:
%      [v1.x v1.y v1.z v2.x v2.y v2.z v3.x v3.y v3.z l1 l2 l3]
%
% v1, v2, v3 are eignvectors and l1, l2, l3 are eigenvalues, in sorted order.
%
% Also see eig2nifti()

loc = strfind(file, '.gz');
len = length(file);

if (~isempty(loc)) && (len-loc(end) == 2)  % .gz extension found
   
   if ~exist('workdir','var')
      workdir = tempname;
      temp_folder = 1;
   end
   niiData = gunzip(file, workdir);
   data = load_untouch_eig(char(niiData));
   delete(char(niiData));
   
   % fix data.fileprefix to reflect the actual folder
   data.fileprefix = remove_extension(file);
   
   % delete temp_dir if created
   if exist('temp_folder','var') && temp_folder == 1 
      rmdir(workdir, 's')
   end
else
   data = load_untouch_eig(file);
end

end


% modified from load_untouch_nii
function nii = load_untouch_eig(filename)

%  Read the dataset header
[nii.hdr,nii.filetype,nii.fileprefix,nii.machine] = load_nii_hdr(filename);

if nii.filetype == 0
   nii.hdr = load_untouch0_nii_hdr(nii.fileprefix,nii.machine);
   nii.ext = [];
else
   nii.hdr = load_untouch_nii_hdr(nii.fileprefix,nii.machine,nii.filetype);
   %  Read the header extension
   nii.ext = load_nii_ext(filename);
end

%  Read the dataset body
[nii.img,nii.hdr] = read_image_eig(nii.hdr, nii.filetype, nii.fileprefix, ...
   nii.machine);

nii.untouch = 1;

return;
end % load_untouch_nii


%---------------------------------------------------------------------
function [img,hdr] = read_image_eig(hdr,filetype,fileprefix,machine)

switch filetype
   case {0, 1}
      error('filetype should be 2.') %fn = [fileprefix '.img'];
   case 2
      fn = [fileprefix '.nii'];
end

fid = fopen(fn,'r',machine);

if fid < 0,
   msg = sprintf('Cannot open file %s.',fn);
   error(msg);
end

% Set for .eig file
%
%  hdr.dime.datatype = int16(0);
%  hdr.dime.bitpix = int16(384);
%  precision = 'float32';

switch hdr.dime.datatype
   case   0,
      hdr.dime.bitpix = 384;
      precision = 'float32';
   otherwise
      error('This file does not looks like a .eig file saved by BDP. \nDatatype is not zero: hdr.dime.datatype = %d', hdr.dime.datatype);
end

% eig file check 
if hdr.dime.dim(1)~=4 || hdr.dime.dim(5)~=12
   hdr.dime
   error('This file does not looks like a .eig file saved by BDP.');
end

% compute size of full volume
hdr.dime.dim(2:end) = abs(hdr.dime.dim(2:end));
d1 = hdr.dime.dim(2);
d2 = hdr.dime.dim(3);
d3 = hdr.dime.dim(4);
d4 = hdr.dime.dim(5);
vox_siz = d1*d2*d3*12;

%  move pointer to the start of image block
switch filetype
   case {0, 1}
      fseek(fid, 0, 'bof');
   case 2
      fseek(fid, hdr.dime.vox_offset, 'bof');
end

img = fread(fid, vox_siz, sprintf('*%s',precision));
fclose(fid);

%  Update the global min and max values
hdr.dime.glmax = double(max(img(:)));
hdr.dime.glmin = double(min(img(:)));

% reshape & permute to reflect correct image matrix for .eig
img = (reshape(img, [12 d1 d2 d3]));
img = permute(img, [2 3 4 1]);

end  % read_image

