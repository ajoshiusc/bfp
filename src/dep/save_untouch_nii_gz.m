function [ savedFileName ] = save_untouch_nii_gz(niiData, fileName, datatype, workdir)
% Saves compressed GNU zip NIfTI-1 files.
% WARNING: This function 'silently' overwrites the destination filename.
%
%   niiData: nii structure
%   fileName: Destination filename, with/without '.nii' or '.nii.gz'.
%             Use [] for no filename. When [] is used, filename is
%             retrived from the nii structure or the file is compressed
%             with the same filename.
%   datatype: (Optional) NiFTI datatype; must be integer. (see save_nii.m
%             for exhaustive list) Few common data types:
%                2    Unsigned char; uint8, NIFTI_TYPE_UINT8 (for BrainSuite masks)
%                4    Signed short; int16; NIFTI_TYPE_INT16
%               16    Floating point; single or float32; NIFTI_TYPE_FLOAT32
%               64    Double precision; double or float64; DT_FLOAT64, NIFTI_TYPE_FLOAT64
%   workdir (optional)
% ________________________________________________________________________________
%
% This function uses functions from following toolbox:
%
% Tools for NIfTI and ANALYZE image
% (Code covered by the BSD License)
% http://www.mathworks.com/matlabcentral/fileexchange/8797
% http://research.baycrest.org/~jimmy/NIfTI/
%
% Copyright (c) 2009, Jimmy Shen
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%   * Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in
%     the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%

if isstruct(niiData)
   
   if (nargin==3) && ischar(datatype)
      workdir = datatype;
      clear datatype;
   end
   
   if isempty(fileName)
      [pathstr, fileName] = fileBaseName(niiData.fileprefix);
   else
      [pathstr, fileName] = fileBaseName(remove_extension(fileName));
   end
   
   if isempty(pathstr)
      pathstr = pwd;
   end
   
   if ~exist('workdir', 'var')
      workdir = tempname;
      mkdir(workdir)
      temp_folder = 1;
   end
   
   if exist('datatype','var')
      niiData.hdr.dime.datatype = datatype;
      niiData.hdr.dime.bitpix = datatype2bitpix(datatype);
   end
   
   save_untouch_nii_SAVEUNTOUCHNIIGZ(niiData, fullfile(workdir, [fileName '.nii']));
   fname = gzip(fullfile(workdir, [fileName '.nii']), pathstr);
   savedFileName = fname{1};
   
   if exist('temp_folder','var') && (temp_folder == 1)
      rmdir(workdir, 's');
   end
   
else
   error('niiData must be a nii structure.')
end

end


function save_untouch_nii_SAVEUNTOUCHNIIGZ(nii, filename)
   
   if ~exist('nii','var') | isempty(nii) | ~isfield(nii,'hdr') | ...
	~isfield(nii,'img') | ~exist('filename','var') | isempty(filename)

      error('Usage: save_untouch_nii(nii, filename)');
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

   v = version;

   %  Check file extension. If .gz, unpack it into temp folder
   %
   if length(filename) > 2 & strcmp(filename(end-2:end), '.gz')

      if ~strcmp(filename(end-6:end), '.img.gz') & ...
	 ~strcmp(filename(end-6:end), '.hdr.gz') & ...
	 ~strcmp(filename(end-6:end), '.nii.gz')

         error('Please check filename.');
      end

      if str2num(v(1:3)) < 7.1 | ~usejava('jvm')
         error('Please use MATLAB 7.1 (with java) and above, or run gunzip outside MATLAB.');
      else
         gzFile = 1;
         filename = filename(1:end-3);
      end
   end

   [p,f] = fileparts(filename);
   fileprefix = fullfile(p, f);

   write_nii_SAVEUNTOUCHNIIGZ(nii, filetype, fileprefix);

   %  gzip output file if requested
   %
   if exist('gzFile', 'var')
      if filetype == 1
         gzip([fileprefix, '.img']);
         delete([fileprefix, '.img']);
         gzip([fileprefix, '.hdr']);
         delete([fileprefix, '.hdr']);
      elseif filetype == 2
         gzip([fileprefix, '.nii']);
         delete([fileprefix, '.nii']);
      end;
   end;

%   %  So earlier versions of SPM can also open it with correct originator
 %  %
  % if filetype == 0
   %   M=[[diag(nii.hdr.dime.pixdim(2:4)) -[nii.hdr.hist.originator(1:3).*nii.hdr.dime.pixdim(2:4)]'];[0 0 0 1]];
    %  save(fileprefix, 'M');
%   elseif filetype == 1
 %     M=[];
  %    save(fileprefix, 'M');
   %end
   
   return					% save_untouch_nii
end

function write_nii_SAVEUNTOUCHNIIGZ(nii, filetype, fileprefix)

   hdr = nii.hdr;

   if isfield(nii,'ext') & ~isempty(nii.ext)
      ext = nii.ext;
      [ext, esize_total] = verify_nii_ext_SAVEUNTOUCHNIIGZ(ext);
   else
      ext = [];
   end

   switch double(hdr.dime.datatype),
   case   1,
      hdr.dime.bitpix = int16(1 ); precision = 'ubit1';
   case   2,
      hdr.dime.bitpix = int16(8 ); precision = 'uint8';
   case   4,
      hdr.dime.bitpix = int16(16); precision = 'int16';
   case   8,
      hdr.dime.bitpix = int16(32); precision = 'int32';
   case  16,
      hdr.dime.bitpix = int16(32); precision = 'float32';
   case  32,
      hdr.dime.bitpix = int16(64); precision = 'float32';
   case  64,
      hdr.dime.bitpix = int16(64); precision = 'float64';
   case 128,
      hdr.dime.bitpix = int16(24); precision = 'uint8';
   case 256 
      hdr.dime.bitpix = int16(8 ); precision = 'int8';
   case 512 
      hdr.dime.bitpix = int16(16); precision = 'uint16';
   case 768 
      hdr.dime.bitpix = int16(32); precision = 'uint32';
   case 1024
      hdr.dime.bitpix = int16(64); precision = 'int64';
   case 1280
      hdr.dime.bitpix = int16(64); precision = 'uint64';
   case 1792,
      hdr.dime.bitpix = int16(128); precision = 'float64';
   otherwise
      error('This datatype is not supported');
   end
   
%   hdr.dime.glmax = round(double(max(nii.img(:))));
 %  hdr.dime.glmin = round(double(min(nii.img(:))));
   
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
      save_untouch_nii_hdr_SAVEUNTOUCHNIIGZ(hdr, fid);

      if ~isempty(ext)
         save_nii_ext_SAVEUNTOUCHNIIGZ(ext, fid);
      end
   elseif filetype == 1
      fid = fopen(sprintf('%s.hdr',fileprefix),'w');
      
      if fid < 0,
         msg = sprintf('Cannot open file %s.hdr.',fileprefix);
         error(msg);
      end
      
      hdr.dime.vox_offset = 0;
      hdr.hist.magic = 'ni1';
      save_untouch_nii_hdr_SAVEUNTOUCHNIIGZ(hdr, fid);

      if ~isempty(ext)
         save_nii_ext_SAVEUNTOUCHNIIGZ(ext, fid);
      end
      
      fclose(fid);
      fid = fopen(sprintf('%s.img',fileprefix),'w');
   else
      error('filetype is zero! Saving to ANALYZE 7.5 is not supported.')
      %  fid = fopen(sprintf('%s.hdr',fileprefix),'w');
      %
      %  if fid < 0,
      %     msg = sprintf('Cannot open file %s.hdr.',fileprefix);
      %     error(msg);
      %   end
      %
      %  save_untouch0_nii_hdr(hdr, fid);
      %
      %  fclose(fid);
      %  fid = fopen(sprintf('%s.img',fileprefix),'w');
   end

   ScanDim = double(hdr.dime.dim(5));		% t
   SliceDim = double(hdr.dime.dim(4));		% z
   RowDim   = double(hdr.dime.dim(3));		% y
   PixelDim = double(hdr.dime.dim(2));		% x
   SliceSz  = double(hdr.dime.pixdim(4));
   RowSz    = double(hdr.dime.pixdim(3));
   PixelSz  = double(hdr.dime.pixdim(2));
   
   x = 1:PixelDim;
   
   if filetype == 2 & isempty(ext)
      skip_bytes = double(hdr.dime.vox_offset) - 348;
   else
      skip_bytes = 0;
   end

   if double(hdr.dime.datatype) == 128

      %  RGB planes are expected to be in the 4th dimension of nii.img
      %
      if(size(nii.img,4)~=3)
         error(['The NII structure does not appear to have 3 RGB color planes in the 4th dimension']);
      end

      nii.img = permute(nii.img, [4 1 2 3 5 6 7 8]);
   end

   %  For complex float32 or complex float64, voxel values
   %  include [real, imag]
   %
   if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
      real_img = real(nii.img(:))';
      nii.img = imag(nii.img(:))';
      nii.img = [real_img; nii.img];
   end

   if skip_bytes
      fwrite(fid, zeros(1,skip_bytes), 'uint8');
   end

   fwrite(fid, nii.img, precision);
%   fwrite(fid, nii.img, precision, skip_bytes);        % error using skip
   fclose(fid);

   return;					% write_nii

end

function save_untouch_nii_hdr_SAVEUNTOUCHNIIGZ(hdr, fid)

   if ~isequal(hdr.hk.sizeof_hdr,348),
      error('hdr.hk.sizeof_hdr must be 348.');
   end

   write_header_SAVEUNTOUCHNIIGZ(hdr, fid);

   return;					% save_nii_hdr
end

%---------------------------------------------------------------------
function write_header_SAVEUNTOUCHNIIGZ(hdr, fid)

        %  Original header structures
	%  struct dsr				/* dsr = hdr */
	%       { 
	%       struct header_key hk;            /*   0 +  40       */
	%       struct image_dimension dime;     /*  40 + 108       */
	%       struct data_history hist;        /* 148 + 200       */
	%       };                               /* total= 348 bytes*/
   
   header_key_SAVEUNTOUCHNIIGZ(fid, hdr.hk);
   image_dimension_SAVEUNTOUCHNIIGZ(fid, hdr.dime);
   data_history_SAVEUNTOUCHNIIGZ(fid, hdr.hist);
   
   %  check the file size is 348 bytes
   %
   fbytes = ftell(fid);
   
   if ~isequal(fbytes,348),
      msg = sprintf('Header size is not 348 bytes.');
      warning(msg);
   end
    
   return;					% write_header
end

%---------------------------------------------------------------------
function header_key_SAVEUNTOUCHNIIGZ(fid, hk)
   
   fseek(fid,0,'bof');

	%  Original header structures    
	%  struct header_key                      /* header key      */ 
	%       {                                /* off + size      */
	%       int sizeof_hdr                   /*  0 +  4         */
	%       char data_type[10];              /*  4 + 10         */
	%       char db_name[18];                /* 14 + 18         */
	%       int extents;                     /* 32 +  4         */
	%       short int session_error;         /* 36 +  2         */
	%       char regular;                    /* 38 +  1         */
	%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
	%       };                               /* total=40 bytes  */
        
   fwrite(fid, hk.sizeof_hdr(1),    'int32');	% must be 348.
    
   % data_type = sprintf('%-10s',hk.data_type);	% ensure it is 10 chars from left
   % fwrite(fid, data_type(1:10), 'uchar');
   pad = zeros(1, 10-length(hk.data_type));
   hk.data_type = [hk.data_type  char(pad)];
   fwrite(fid, hk.data_type(1:10), 'uchar');
    
   % db_name   = sprintf('%-18s', hk.db_name);	% ensure it is 18 chars from left
   % fwrite(fid, db_name(1:18), 'uchar');
   pad = zeros(1, 18-length(hk.db_name));
   hk.db_name = [hk.db_name  char(pad)];
   fwrite(fid, hk.db_name(1:18), 'uchar');
    
   fwrite(fid, hk.extents(1),       'int32');
   fwrite(fid, hk.session_error(1), 'int16');
   fwrite(fid, hk.regular(1),       'uchar');	% might be uint8
    
   % fwrite(fid, hk.hkey_un0(1),    'uchar');
   % fwrite(fid, hk.hkey_un0(1),    'uint8');
   fwrite(fid, hk.dim_info(1),      'uchar');
    
   return;					% header_key
end

%---------------------------------------------------------------------
function image_dimension_SAVEUNTOUCHNIIGZ(fid, dime)

	%  Original header structures        
	%  struct image_dimension
	%       {                                /* off + size      */
	%       short int dim[8];                /* 0 + 16          */
	%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
	%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
	%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
	%       short int intent_code;   % short int unused1;   /* 28 + 2 */
	%       short int datatype;              /* 30 + 2          */
	%       short int bitpix;                /* 32 + 2          */
	%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
	%       float pixdim[8];                 /* 36 + 32         */
	%			/*
	%				pixdim[] specifies the voxel dimensions:
	%				pixdim[1] - voxel width
	%				pixdim[2] - voxel height
	%				pixdim[3] - interslice distance
	%				pixdim[4] - volume timing, in msec
	%					..etc
	%			*/
	%       float vox_offset;                /* 68 + 4          */
	%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
	%       float scl_inter;   % float funused1;      /* 76 + 4 */
	%       short slice_end;   % float funused2;      /* 80 + 2 */
	%       char slice_code;   % float funused2;      /* 82 + 1 */
	%       char xyzt_units;   % float funused2;      /* 83 + 1 */
	%       float cal_max;                   /* 84 + 4          */
	%       float cal_min;                   /* 88 + 4          */
	%       float slice_duration;   % int compressed; /* 92 + 4 */
	%       float toffset;   % int verified;          /* 96 + 4 */
	%       int glmax;                       /* 100 + 4         */
	%       int glmin;                       /* 104 + 4         */
	%       };                               /* total=108 bytes */
	
   fwrite(fid, dime.dim(1:8),        'int16');
   fwrite(fid, dime.intent_p1(1),  'float32');
   fwrite(fid, dime.intent_p2(1),  'float32');
   fwrite(fid, dime.intent_p3(1),  'float32');
   fwrite(fid, dime.intent_code(1),  'int16');
   fwrite(fid, dime.datatype(1),     'int16');
   fwrite(fid, dime.bitpix(1),       'int16');
   fwrite(fid, dime.slice_start(1),  'int16');
   fwrite(fid, dime.pixdim(1:8),   'float32');
   fwrite(fid, dime.vox_offset(1), 'float32');
   fwrite(fid, dime.scl_slope(1),  'float32');
   fwrite(fid, dime.scl_inter(1),  'float32');
   fwrite(fid, dime.slice_end(1),    'int16');
   fwrite(fid, dime.slice_code(1),   'uchar');
   fwrite(fid, dime.xyzt_units(1),   'uchar');
   fwrite(fid, dime.cal_max(1),    'float32');
   fwrite(fid, dime.cal_min(1),    'float32');
   fwrite(fid, dime.slice_duration(1), 'float32');
   fwrite(fid, dime.toffset(1),    'float32');
   fwrite(fid, dime.glmax(1),        'int32');
   fwrite(fid, dime.glmin(1),        'int32');
   
   return;					% image_dimension
end

%---------------------------------------------------------------------
function data_history_SAVEUNTOUCHNIIGZ(fid, hist)
    
	% Original header structures
	%struct data_history       
	%       {                                /* off + size      */
	%       char descrip[80];                /* 0 + 80          */
	%       char aux_file[24];               /* 80 + 24         */
	%       short int qform_code;            /* 104 + 2         */
	%       short int sform_code;            /* 106 + 2         */
	%       float quatern_b;                 /* 108 + 4         */
	%       float quatern_c;                 /* 112 + 4         */
	%       float quatern_d;                 /* 116 + 4         */
	%       float qoffset_x;                 /* 120 + 4         */
	%       float qoffset_y;                 /* 124 + 4         */
	%       float qoffset_z;                 /* 128 + 4         */
	%       float srow_x[4];                 /* 132 + 16        */
	%       float srow_y[4];                 /* 148 + 16        */
	%       float srow_z[4];                 /* 164 + 16        */
	%       char intent_name[16];            /* 180 + 16        */
	%       char magic[4];   % int smin;     /* 196 + 4         */
	%       };                               /* total=200 bytes */
	
   % descrip     = sprintf('%-80s', hist.descrip);     % 80 chars from left
   % fwrite(fid, descrip(1:80),    'uchar');
   pad = zeros(1, 80-length(hist.descrip));
   hist.descrip = [hist.descrip  char(pad)];
   fwrite(fid, hist.descrip(1:80), 'uchar');
    
   % aux_file    = sprintf('%-24s', hist.aux_file);    % 24 chars from left
   % fwrite(fid, aux_file(1:24),   'uchar');
   pad = zeros(1, 24-length(hist.aux_file));
   hist.aux_file = [hist.aux_file  char(pad)];
   fwrite(fid, hist.aux_file(1:24), 'uchar');
    
   fwrite(fid, hist.qform_code,    'int16');
   fwrite(fid, hist.sform_code,    'int16');
   fwrite(fid, hist.quatern_b,   'float32');
   fwrite(fid, hist.quatern_c,   'float32');
   fwrite(fid, hist.quatern_d,   'float32');
   fwrite(fid, hist.qoffset_x,   'float32');
   fwrite(fid, hist.qoffset_y,   'float32');
   fwrite(fid, hist.qoffset_z,   'float32');
   fwrite(fid, hist.srow_x(1:4), 'float32');
   fwrite(fid, hist.srow_y(1:4), 'float32');
   fwrite(fid, hist.srow_z(1:4), 'float32');

   % intent_name = sprintf('%-16s', hist.intent_name);	% 16 chars from left
   % fwrite(fid, intent_name(1:16),    'uchar');
   pad = zeros(1, 16-length(hist.intent_name));
   hist.intent_name = [hist.intent_name  char(pad)];
   fwrite(fid, hist.intent_name(1:16), 'uchar');
    
   % magic	= sprintf('%-4s', hist.magic);		% 4 chars from left
   % fwrite(fid, magic(1:4),           'uchar');
   pad = zeros(1, 4-length(hist.magic));
   hist.magic = [hist.magic  char(pad)];
   fwrite(fid, hist.magic(1:4),        'uchar');
    
   return;					% data_history
end

function save_nii_ext_SAVEUNTOUCHNIIGZ(ext, fid)

   if ~exist('ext','var') | ~exist('fid','var')
      error('Usage: save_nii_ext(ext, fid)');
   end

   if ~isfield(ext,'extension') | ~isfield(ext,'section') | ~isfield(ext,'num_ext')
      error('Wrong header extension');
   end

   write_ext_SAVEUNTOUCHNIIGZ(ext, fid);

   return;                                      % save_nii_ext
end

%---------------------------------------------------------------------
function write_ext_SAVEUNTOUCHNIIGZ(ext, fid)

   fwrite(fid, ext.extension, 'uchar');

   for i=1:ext.num_ext
      fwrite(fid, ext.section(i).esize, 'int32');
      fwrite(fid, ext.section(i).ecode, 'int32');
      fwrite(fid, ext.section(i).edata, 'uchar');
   end

   return;                                      % write_ext
end

function [ext, esize_total] = verify_nii_ext_SAVEUNTOUCHNIIGZ(ext)

   if ~isfield(ext, 'section')
      error('Incorrect NIFTI header extension structure.');
   elseif ~isfield(ext, 'num_ext')
      ext.num_ext = length(ext.section);
   elseif ~isfield(ext, 'extension')
      ext.extension = [1 0 0 0];
   end

   esize_total = 0;

   for i=1:ext.num_ext
      if ~isfield(ext.section(i), 'ecode') | ~isfield(ext.section(i), 'edata')
         error('Incorrect NIFTI header extension structure.');
      end

      ext.section(i).esize = ceil((length(ext.section(i).edata)+8)/16)*16;
      ext.section(i).edata = ...
	[ext.section(i).edata ...
	 zeros(1,ext.section(i).esize-length(ext.section(i).edata)-8)];
      esize_total = esize_total + ext.section(i).esize;
   end

   return                                       % verify_nii_ext
end
