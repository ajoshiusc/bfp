function generate_eig_file(L1, L2, L3, V1, V2, V3, outputFileName)
%GENERATE_EIG Generate .eig file. Eigenvalues (L1, L2, L3) &
% Eigenvectors(V1, V2, V3) can be nii structures or .nii.gz or .nii files.
%
%   outputFileName - file name with/without extension (with full path)
%

workdir = tempname();
mkdir(workdir)

if ischar(V1)
   V1 = load_untouch_nii_gz(V1, true, workdir);
elseif ~isfield(V1,'untouch') || V1.untouch ~= 1
   error('Please use ''load_untouch_nii.m'' to load the structure.');
end

num_vox = prod(V1.hdr.dime.dim(2:4));
eig = V1;
eig.img = [];
eig.hdr.dime.dim(5) = 12;
eig.hdr.hk.data_type = 0;
eig.hdr.dime.bitpix = 384;
eig.hdr.hist.descrip = '.eig file for BrainSuite';


eig.img = single(reshape(V1.img(:), [num_vox 3]));
clear V1


if ischar(V2)
   V2 = load_untouch_nii_gz(V2, true, workdir);
elseif ~isfield(V2,'untouch') || V2.untouch ~= 1
   error('Please use ''load_untouch_nii.m'' to load the structure.');
end
eig.img = [eig.img single(reshape(V2.img(:), [num_vox 3]))];
clear V2


if ischar(V3)
   V3 = load_untouch_nii_gz(V3, true, workdir);
elseif ~isfield(V3,'untouch') || V3.untouch ~= 1
   error('Please use ''load_untouch_nii.m'' to load the structure.');
end
eig.img = [eig.img single(reshape(V3.img(:), [num_vox 3]))];
clear V3


if ischar(L1)
   L1 = load_untouch_nii_gz(L1, true, workdir);
elseif ~isfield(L1,'untouch') || L1.untouch ~= 1
   error('Please use ''load_untouch_nii.m'' to load the structure.');
end
eig.img = [eig.img single(L1.img(:))];
clear L1


if ischar(L2)
   L2 = load_untouch_nii_gz(L2, true, workdir);
elseif ~isfield(L2,'untouch') || L2.untouch ~= 1
   error('Please use ''load_untouch_nii.m'' to load the structure.');
end
eig.img = [eig.img single(L2.img(:))];
clear L2


if ischar(L3)
   L3 = load_untouch_nii_gz(L3, true, workdir);
elseif ~isfield(L3,'untouch') || L3.untouch ~= 1
   error('Please use ''load_untouch_nii.m'' to load the structure.');
end
eig.img = [eig.img single(L3.img(:))];
clear L3


% correct order of eigen-data - 12*num_vox
eig.img = eig.img';
eig.img = eig.img(:);

[pathstr, file_name] = fileBaseName(remove_extension(outputFileName));
if isempty(pathstr)
   pathstr = pwd;
end

tempfile = [workdir '/' file_name '.eig.nii'];
save_untouch_eig(eig, tempfile);
gzip(tempfile, pathstr);


rmdir(workdir, 's');
end



% modified from save_untouch_nii
function save_untouch_eig(nii, filename)

if ~exist('nii','var') | isempty(nii) | ~isfield(nii,'hdr') | ...
      ~isfield(nii,'img') | ~exist('filename','var') | isempty(filename)
   
   error('nii structure does not look correct. Check all nii structure.');
end

if ~isfield(nii,'untouch') | nii.untouch == 0
   error('Usage: please use ''save_nii.m'' for the modified structure.');
end

if isfield(nii.hdr.hist,'magic') & strcmp(nii.hdr.hist.magic(1:3),'ni1')
   filetype = 1;
elseif isfield(nii.hdr.hist,'magic') & strcmp(nii.hdr.hist.magic(1:3),'n+1')
   filetype = 2;
else
   filetype = 0;
end

[p,f] = fileparts(filename);
fileprefix = fullfile(p, f);

write_nii_SAVE_UNTOUCH_EIG(nii, filetype, fileprefix);

end


function write_nii_SAVE_UNTOUCH_EIG(nii, filetype, fileprefix)

hdr = nii.hdr;

if isfield(nii,'ext') && ~isempty(nii.ext)
   ext = nii.ext;
   [ext, esize_total] = verify_nii_ext(ext);
else
   ext = [];
end

% set for eig files
hdr.dime.datatype = int16(0);
hdr.dime.bitpix = int16(384);
precision = 'float32';

if filetype == 2
   fid = fopen(sprintf('%s.nii',fileprefix),'w');
   
   if fid < 0,
      msg = sprintf('Cannot open file %s.nii.',fileprefix);
      error(msg);
   end
   
   hdr.dime.vox_offset = 352;
   
   if ~isempty(ext)
      hdr.dime.vox_offset = hdr.dime.vox_offset + esize_total;
   end
   
   hdr.hist.magic = 'n+1';
   save_untouch_nii_hdr(hdr, fid);
   
   if ~isempty(ext)
      save_nii_ext(ext, fid);
   end
elseif filetype == 1
   error ('filetype should not be 1. Analyze format is not supported.')
else
   error ('filetype should be 2. Any other format is not supported.')
end

if filetype == 2 && isempty(ext)
   skip_bytes = double(hdr.dime.vox_offset) - 348;
else
   skip_bytes = 0;
end


if skip_bytes
   fwrite(fid, zeros(1,skip_bytes), 'uint8');
end

fwrite(fid, nii.img, precision);
fclose(fid);

end



