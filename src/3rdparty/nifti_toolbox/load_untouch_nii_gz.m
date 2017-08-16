function [data, orig_scl_slope, orig_scl_inter, orig_datatype] = load_untouch_nii_gz(file, scl_flag, workdir)
% Wrapper around load_untouch_nii_hdr to work with .gz files. workdir, when specified, is the
% temporary directory to be used for unzipping. When workdir is not specified then a
% random-temporary folder is used. When scl_flag is true the img data is scaled appropriately, if
% required by the header, and hdr parameters are adjusted accordingly (orig_scl_slope,
% orig_scl_inter, orig_datatype can be used to access the original values of these parameters). By
% default, no scaling is performed.
%
% Usage: 
%     [data] = load_untouch_nii_gz(file)
%     [data] = load_untouch_nii_gz(file, scl_flag)
%     [data] = load_untouch_nii_gz(file, workdir)
%     [data] = load_untouch_nii_gz(file, scl_flag, workdir)
%     [data, orig_scl_slope, orig_scl_inter, orig_datatype] = load_untouch_nii_gz(...)
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

if nargin==1
   workdir = tempname();
   temp_folder = 1;
   scl_flag = false;
   
elseif nargin==2 && islogical(scl_flag)
   workdir = tempname();
   temp_folder = 1;   
   
elseif nargin==2 && ~islogical(scl_flag)
   workdir = scl_flag;
   scl_flag = false;
   temp_folder = 0;
else
   temp_folder = 0;
end

loc = strfind(file, '.gz');
len = length(file);

if (~isempty(loc)) && (len-loc(end) == 2)  % .gz extension found

   niiData = gunzip(file, workdir);
   data = load_untouch_nii_LOADUNTOUCHNIIGZ(char(niiData));
   delete(char(niiData));
   
   % fix data.fileprefix to reflect the actual folder
   data.fileprefix = remove_extension(file);
   
   if temp_folder == 1 % delete temp_dir if created
      rmdir(workdir, 's')
   end
   
else
   data = load_untouch_nii_LOADUNTOUCHNIIGZ(file);
end

% scale data if scl_slope is non-zero
orig_scl_slope = double(data.hdr.dime.scl_slope);
orig_scl_inter = double(data.hdr.dime.scl_inter);
orig_datatype = data.hdr.dime.datatype;
if scl_flag && (data.hdr.dime.scl_slope ~= 0)
   if (data.hdr.dime.scl_slope == 1.0 && data.hdr.dime.scl_inter == 0 )
      data.hdr.dime.scl_slope = 0;
   else
      data.img = (double(data.img) * double(data.hdr.dime.scl_slope)) + double(data.hdr.dime.scl_inter);
      data.hdr.dime.scl_slope = 0;
      data.hdr.dime.scl_inter = 0;
      data.hdr.dime.datatype = 16; % single or float32
   end
end

end



function nii = load_untouch_nii_LOADUNTOUCHNIIGZ(filename, img_idx, dim5_idx, dim6_idx, dim7_idx, ...
   old_RGB, slice_idx)

if ~exist('filename','var')
   error('Usage: nii = load_untouch_nii(filename, [img_idx], [dim5_idx], [dim6_idx], [dim7_idx], [old_RGB], [slice_idx])');
end

if ~exist('img_idx','var') | isempty(img_idx)
   img_idx = [];
end

if ~exist('dim5_idx','var') | isempty(dim5_idx)
   dim5_idx = [];
end

if ~exist('dim6_idx','var') | isempty(dim6_idx)
   dim6_idx = [];
end

if ~exist('dim7_idx','var') | isempty(dim7_idx)
   dim7_idx = [];
end

if ~exist('old_RGB','var') | isempty(old_RGB)
   old_RGB = 0;
end

if ~exist('slice_idx','var') | isempty(slice_idx)
   slice_idx = [];
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
   elseif strcmp(filename(end-6:end), '.img.gz')
      filename1 = filename;
      filename2 = filename;
      filename2(end-6:end) = '';
      filename2 = [filename2, '.hdr.gz'];
      
      tmpDir = tempname;
      mkdir(tmpDir);
      gzFileName = filename;
      
      filename1 = gunzip(filename1, tmpDir);
      filename2 = gunzip(filename2, tmpDir);
      filename = char(filename1);	% convert from cell to string
   elseif strcmp(filename(end-6:end), '.hdr.gz')
      filename1 = filename;
      filename2 = filename;
      filename2(end-6:end) = '';
      filename2 = [filename2, '.img.gz'];
      
      tmpDir = tempname;
      mkdir(tmpDir);
      gzFileName = filename;
      
      filename1 = gunzip(filename1, tmpDir);
      filename2 = gunzip(filename2, tmpDir);
      filename = char(filename1);	% convert from cell to string
   elseif strcmp(filename(end-6:end), '.nii.gz')
      tmpDir = tempname;
      mkdir(tmpDir);
      gzFileName = filename;
      filename = gunzip(filename, tmpDir);
      filename = char(filename);	% convert from cell to string
   end
end

%  Read the dataset header
%
[nii.hdr,nii.filetype,nii.fileprefix,nii.machine] = load_nii_hdr_LOADUNTOUCHNIIGZ(filename);

if nii.filetype == 0
   nii.hdr = load_untouch0_nii_hdr_LOADUNTOUCHNIIGZ(nii.fileprefix,nii.machine);
   nii.ext = [];
else
   nii.hdr = load_untouch_nii_hdr_LOADUNTOUCHNIIGZ(nii.fileprefix,nii.machine,nii.filetype);
   
   %  Read the header extension
   %
   nii.ext = load_nii_ext_LOADUNTOUCHNIIGZ(filename);
end

%  Read the dataset body
%
[nii.img,nii.hdr] = load_untouch_nii_img_LOADUNTOUCHNIIGZ(nii.hdr,nii.filetype,nii.fileprefix, ...
   nii.machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB,slice_idx);

%  Perform some of sform/qform transform
%
%   nii = xform_nii(nii, tolerance, preferredForm);

nii.untouch = 1;


%  Clean up after gunzip
%
if exist('gzFileName', 'var')
   
   %  fix fileprefix so it doesn't point to temp location
   %
   nii.fileprefix = gzFileName(1:end-7);
   rmdir(tmpDir,'s');
end


return					% load_untouch_nii
end


function ext = load_nii_ext_LOADUNTOUCHNIIGZ(filename)

if ~exist('filename','var'),
   error('Usage: ext = load_nii_ext(filename)');
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
   elseif strcmp(filename(end-6:end), '.img.gz')
      filename1 = filename;
      filename2 = filename;
      filename2(end-6:end) = '';
      filename2 = [filename2, '.hdr.gz'];
      
      tmpDir = tempname;
      mkdir(tmpDir);
      gzFileName = filename;
      
      filename1 = gunzip(filename1, tmpDir);
      filename2 = gunzip(filename2, tmpDir);
      filename = char(filename1);	% convert from cell to string
   elseif strcmp(filename(end-6:end), '.hdr.gz')
      filename1 = filename;
      filename2 = filename;
      filename2(end-6:end) = '';
      filename2 = [filename2, '.img.gz'];
      
      tmpDir = tempname;
      mkdir(tmpDir);
      gzFileName = filename;
      
      filename1 = gunzip(filename1, tmpDir);
      filename2 = gunzip(filename2, tmpDir);
      filename = char(filename1);	% convert from cell to string
   elseif strcmp(filename(end-6:end), '.nii.gz')
      tmpDir = tempname;
      mkdir(tmpDir);
      gzFileName = filename;
      filename = gunzip(filename, tmpDir);
      filename = char(filename);	% convert from cell to string
   end
end

machine = 'ieee-le';
new_ext = 0;

if findstr('.nii',filename) & strcmp(filename(end-3:end), '.nii')
   new_ext = 1;
   filename(end-3:end)='';
end

if findstr('.hdr',filename) & strcmp(filename(end-3:end), '.hdr')
   filename(end-3:end)='';
end

if findstr('.img',filename) & strcmp(filename(end-3:end), '.img')
   filename(end-3:end)='';
end

if new_ext
   fn = sprintf('%s.nii',filename);
   
   if ~exist(fn)
      msg = sprintf('Cannot find file "%s.nii".', filename);
      error(msg);
   end
else
   fn = sprintf('%s.hdr',filename);
   
   if ~exist(fn)
      msg = sprintf('Cannot find file "%s.hdr".', filename);
      error(msg);
   end
end

fid = fopen(fn,'r',machine);
vox_offset = 0;

if fid < 0,
   msg = sprintf('Cannot open file %s.',fn);
   error(msg);
else
   fseek(fid,0,'bof');
   
   if fread(fid,1,'int32') == 348
      if new_ext
         fseek(fid,108,'bof');
         vox_offset = fread(fid,1,'float32');
      end
      
      ext = read_extension_LOADUNTOUCHNIIGZ(fid, vox_offset);
      fclose(fid);
   else
      fclose(fid);
      
      %  first try reading the opposite endian to 'machine'
      %
      switch machine,
         case 'ieee-le', machine = 'ieee-be';
         case 'ieee-be', machine = 'ieee-le';
      end
      
      fid = fopen(fn,'r',machine);
      
      if fid < 0,
         msg = sprintf('Cannot open file %s.',fn);
         error(msg);
      else
         fseek(fid,0,'bof');
         
         if fread(fid,1,'int32') ~= 348
            
            %  Now throw an error
            %
            msg = sprintf('File "%s" is corrupted.',fn);
            error(msg);
         end
         
         if new_ext
            fseek(fid,108,'bof');
            vox_offset = fread(fid,1,'float32');
         end
         
         ext = read_extension_LOADUNTOUCHNIIGZ(fid, vox_offset);
         fclose(fid);
      end
   end
end


%  Clean up after gunzip
%
if exist('gzFileName', 'var')
   rmdir(tmpDir,'s');
end


return                                       % load_nii_ext
end

%---------------------------------------------------------------------
function ext = read_extension_LOADUNTOUCHNIIGZ(fid, vox_offset)

ext = [];

if vox_offset
   end_of_ext = vox_offset;
else
   fseek(fid, 0, 'eof');
   end_of_ext = ftell(fid);
end

if end_of_ext > 352
   fseek(fid, 348, 'bof');
   ext.extension = fread(fid,4)';
end

if isempty(ext) | ext.extension(1) == 0
   ext = [];
   return;
end

i = 1;

while(ftell(fid) < end_of_ext)
   ext.section(i).esize = fread(fid,1,'int32');
   ext.section(i).ecode = fread(fid,1,'int32');
   ext.section(i).edata = char(fread(fid,ext.section(i).esize-8)');
   i = i + 1;
end

ext.num_ext = length(ext.section);

return                                               % read_extension
end


function [hdr, filetype, fileprefix, machine] = load_nii_hdr_LOADUNTOUCHNIIGZ(fileprefix)

if ~exist('fileprefix','var'),
   error('Usage: [hdr, filetype, fileprefix, machine] = load_nii_hdr(filename)');
end

machine = 'ieee-le';
new_ext = 0;

if findstr('.nii',fileprefix) & strcmp(fileprefix(end-3:end), '.nii')
   new_ext = 1;
   fileprefix(end-3:end)='';
end

if findstr('.hdr',fileprefix) & strcmp(fileprefix(end-3:end), '.hdr')
   fileprefix(end-3:end)='';
end

if findstr('.img',fileprefix) & strcmp(fileprefix(end-3:end), '.img')
   fileprefix(end-3:end)='';
end

if new_ext
   fn = sprintf('%s.nii',fileprefix);
   
   if ~exist(fn)
      msg = sprintf('Cannot find file "%s.nii".', fileprefix);
      error(msg);
   end
else
   fn = sprintf('%s.hdr',fileprefix);
   
   if ~exist(fn)
      msg = sprintf('Cannot find file "%s.hdr".', fileprefix);
      error(msg);
   end
end

fid = fopen(fn,'r',machine);

if fid < 0,
   msg = sprintf('Cannot open file %s.',fn);
   error(msg);
else
   fseek(fid,0,'bof');
   
   if fread(fid,1,'int32') == 348
      hdr = read_header_LOADUNTOUCHNIIGZ(fid);
      fclose(fid);
   else
      fclose(fid);
      
      %  first try reading the opposite endian to 'machine'
      %
      switch machine,
         case 'ieee-le', machine = 'ieee-be';
         case 'ieee-be', machine = 'ieee-le';
      end
      
      fid = fopen(fn,'r',machine);
      
      if fid < 0,
         msg = sprintf('Cannot open file %s.',fn);
         error(msg);
      else
         fseek(fid,0,'bof');
         
         if fread(fid,1,'int32') ~= 348
            
            %  Now throw an error
            %
            msg = sprintf('File "%s" is corrupted.',fn);
            error(msg);
         end
         
         hdr = read_header_LOADUNTOUCHNIIGZ(fid);
         fclose(fid);
      end
   end
end

if strcmp(hdr.hist.magic, 'n+1')
   filetype = 2;
elseif strcmp(hdr.hist.magic, 'ni1')
   filetype = 1;
else
   filetype = 0;
end

return					% load_nii_hdr
end

%---------------------------------------------------------------------
function [ dsr ] = read_header_LOADUNTOUCHNIIGZ(fid)

%  Original header structures
%  struct dsr
%       {
%       struct header_key hk;            /*   0 +  40       */
%       struct image_dimension dime;     /*  40 + 108       */
%       struct data_history hist;        /* 148 + 200       */
%       };                               /* total= 348 bytes*/

dsr.hk   = header_key_LOADUNTOUCHNIIGZ(fid);
dsr.dime = image_dimension_LOADUNTOUCHNIIGZ(fid);
dsr.hist = data_history_LOADUNTOUCHNIIGZ(fid);

%  For Analyze data format
%
if ~strcmp(dsr.hist.magic, 'n+1') & ~strcmp(dsr.hist.magic, 'ni1')
   dsr.hist.qform_code = 0;
   dsr.hist.sform_code = 0;
end

return					% read_header
end

%---------------------------------------------------------------------
function [ hk ] = header_key_LOADUNTOUCHNIIGZ(fid)

fseek(fid,0,'bof');

%  Original header structures
%  struct header_key                     /* header key      */
%       {                                /* off + size      */
%       int sizeof_hdr                   /*  0 +  4         */
%       char data_type[10];              /*  4 + 10         */
%       char db_name[18];                /* 14 + 18         */
%       int extents;                     /* 32 +  4         */
%       short int session_error;         /* 36 +  2         */
%       char regular;                    /* 38 +  1         */
%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
%       };                               /* total=40 bytes  */
%
% int sizeof_header   Should be 348.
% char regular        Must be 'r' to indicate that all images and
%                     volumes are the same size.

v6 = version;
if str2num(v6(1))<6
   directchar = '*char';
else
   directchar = 'uchar=>char';
end

hk.sizeof_hdr    = fread(fid, 1,'int32')';	% should be 348!
hk.data_type     = deblank(fread(fid,10,directchar)');
hk.db_name       = deblank(fread(fid,18,directchar)');
hk.extents       = fread(fid, 1,'int32')';
hk.session_error = fread(fid, 1,'int16')';
hk.regular       = fread(fid, 1,directchar)';
hk.dim_info      = fread(fid, 1,'uchar')';

return					% header_key
end

%---------------------------------------------------------------------
function [ dime ] = image_dimension_LOADUNTOUCHNIIGZ(fid)

%  Original header structures
%  struct image_dimension
%       {                                /* off + size      */
%       short int dim[8];                /* 0 + 16          */
%       /*
%           dim[0]      Number of dimensions in database; usually 4.
%           dim[1]      Image X dimension;  number of *pixels* in an image row.
%           dim[2]      Image Y dimension;  number of *pixel rows* in slice.
%           dim[3]      Volume Z dimension; number of *slices* in a volume.
%           dim[4]      Time points; number of volumes in database
%       */
%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
%       short int intent_code;   % short int unused1;   /* 28 + 2 */
%       short int datatype;              /* 30 + 2          */
%       short int bitpix;                /* 32 + 2          */
%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
%       float pixdim[8];                 /* 36 + 32         */
%	/*
%		pixdim[] specifies the voxel dimensions:
%		pixdim[1] - voxel width, mm
%		pixdim[2] - voxel height, mm
%		pixdim[3] - slice thickness, mm
%		pixdim[4] - volume timing, in msec
%					..etc
%	*/
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

dime.dim        = fread(fid,8,'int16')';
dime.intent_p1  = fread(fid,1,'float32')';
dime.intent_p2  = fread(fid,1,'float32')';
dime.intent_p3  = fread(fid,1,'float32')';
dime.intent_code = fread(fid,1,'int16')';
dime.datatype   = fread(fid,1,'int16')';
dime.bitpix     = fread(fid,1,'int16')';
dime.slice_start = fread(fid,1,'int16')';
dime.pixdim     = fread(fid,8,'float32')';
dime.vox_offset = fread(fid,1,'float32')';
dime.scl_slope  = fread(fid,1,'float32')';
dime.scl_inter  = fread(fid,1,'float32')';
dime.slice_end  = fread(fid,1,'int16')';
dime.slice_code = fread(fid,1,'uchar')';
dime.xyzt_units = fread(fid,1,'uchar')';
dime.cal_max    = fread(fid,1,'float32')';
dime.cal_min    = fread(fid,1,'float32')';
dime.slice_duration = fread(fid,1,'float32')';
dime.toffset    = fread(fid,1,'float32')';
dime.glmax      = fread(fid,1,'int32')';
dime.glmin      = fread(fid,1,'int32')';

return					% image_dimension
end

%---------------------------------------------------------------------
function [ hist ] = data_history_LOADUNTOUCHNIIGZ(fid)

%  Original header structures
%  struct data_history
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

v6 = version;
if str2num(v6(1))<6
   directchar = '*char';
else
   directchar = 'uchar=>char';
end

hist.descrip     = deblank(fread(fid,80,directchar)');
hist.aux_file    = deblank(fread(fid,24,directchar)');
hist.qform_code  = fread(fid,1,'int16')';
hist.sform_code  = fread(fid,1,'int16')';
hist.quatern_b   = fread(fid,1,'float32')';
hist.quatern_c   = fread(fid,1,'float32')';
hist.quatern_d   = fread(fid,1,'float32')';
hist.qoffset_x   = fread(fid,1,'float32')';
hist.qoffset_y   = fread(fid,1,'float32')';
hist.qoffset_z   = fread(fid,1,'float32')';
hist.srow_x      = fread(fid,4,'float32')';
hist.srow_y      = fread(fid,4,'float32')';
hist.srow_z      = fread(fid,4,'float32')';
hist.intent_name = deblank(fread(fid,16,directchar)');
hist.magic       = deblank(fread(fid,4,directchar)');

fseek(fid,253,'bof');
hist.originator  = fread(fid, 5,'int16')';

return					% data_history
end


function hdr = load_untouch0_nii_hdr_LOADUNTOUCHNIIGZ(fileprefix, machine)

fn = sprintf('%s.hdr',fileprefix);
fid = fopen(fn,'r',machine);

if fid < 0,
   msg = sprintf('Cannot open file %s.',fn);
   error(msg);
else
   fseek(fid,0,'bof');
   hdr = read_header_untouch0_LOADUNTOUCHNIIGZ(fid);
   fclose(fid);
end

return					% load_untouch0_nii_hdr
end

%---------------------------------------------------------------------
function [ dsr ] = read_header_untouch0_LOADUNTOUCHNIIGZ(fid)

%  Original header structures
%  struct dsr
%       {
%       struct header_key hk;            /*   0 +  40       */
%       struct image_dimension dime;     /*  40 + 108       */
%       struct data_history hist;        /* 148 + 200       */
%       };                               /* total= 348 bytes*/

dsr.hk   = header_key_untouch0_LOADUNTOUCHNIIGZ(fid);
dsr.dime = image_dimension_untouch0_LOADUNTOUCHNIIGZ(fid);
dsr.hist = data_history_untouch0_LOADUNTOUCHNIIGZ(fid);

return					% read_header
end

%---------------------------------------------------------------------
function [ hk ] = header_key_untouch0_LOADUNTOUCHNIIGZ(fid)

fseek(fid,0,'bof');

%  Original header structures
%  struct header_key                     /* header key      */
%       {                                /* off + size      */
%       int sizeof_hdr                   /*  0 +  4         */
%       char data_type[10];              /*  4 + 10         */
%       char db_name[18];                /* 14 + 18         */
%       int extents;                     /* 32 +  4         */
%       short int session_error;         /* 36 +  2         */
%       char regular;                    /* 38 +  1         */
%       char hkey_un0;                   /* 39 +  1 */
%       };                               /* total=40 bytes  */
%
% int sizeof_header   Should be 348.
% char regular        Must be 'r' to indicate that all images and
%                     volumes are the same size.

v6 = version;
if str2num(v6(1))<6
   directchar = '*char';
else
   directchar = 'uchar=>char';
end

hk.sizeof_hdr    = fread(fid, 1,'int32')';	% should be 348!
hk.data_type     = deblank(fread(fid,10,directchar)');
hk.db_name       = deblank(fread(fid,18,directchar)');
hk.extents       = fread(fid, 1,'int32')';
hk.session_error = fread(fid, 1,'int16')';
hk.regular       = fread(fid, 1,directchar)';
hk.hkey_un0      = fread(fid, 1,directchar)';

return					% header_key
end

%---------------------------------------------------------------------
function [ dime ] = image_dimension_untouch0_LOADUNTOUCHNIIGZ(fid)

%struct image_dimension
%       {                                /* off + size      */
%       short int dim[8];                /* 0 + 16          */
%           /*
%           dim[0]      Number of dimensions in database; usually 4.
%           dim[1]      Image X dimension;  number of *pixels* in an image row.
%           dim[2]      Image Y dimension;  number of *pixel rows* in slice.
%           dim[3]      Volume Z dimension; number of *slices* in a volume.
%           dim[4]      Time points; number of volumes in database
%           */
%       char vox_units[4];               /* 16 + 4          */
%       char cal_units[8];               /* 20 + 8          */
%       short int unused1;               /* 28 + 2          */
%       short int datatype;              /* 30 + 2          */
%       short int bitpix;                /* 32 + 2          */
%       short int dim_un0;               /* 34 + 2          */
%       float pixdim[8];                 /* 36 + 32         */
%			/*
%				pixdim[] specifies the voxel dimensions:
%				pixdim[1] - voxel width, mm
%				pixdim[2] - voxel height, mm
%				pixdim[3] - slice thickness, mm
%               pixdim[4] - volume timing, in msec
%					..etc
%			*/
%       float vox_offset;                /* 68 + 4          */
%       float roi_scale;                 /* 72 + 4          */
%       float funused1;                  /* 76 + 4          */
%       float funused2;                  /* 80 + 4          */
%       float cal_max;                   /* 84 + 4          */
%       float cal_min;                   /* 88 + 4          */
%       int compressed;                  /* 92 + 4          */
%       int verified;                    /* 96 + 4          */
%       int glmax;                       /* 100 + 4         */
%       int glmin;                       /* 104 + 4         */
%       };                               /* total=108 bytes */

v6 = version;
if str2num(v6(1))<6
   directchar = '*char';
else
   directchar = 'uchar=>char';
end

dime.dim        = fread(fid,8,'int16')';
dime.vox_units  = deblank(fread(fid,4,directchar)');
dime.cal_units  = deblank(fread(fid,8,directchar)');
dime.unused1    = fread(fid,1,'int16')';
dime.datatype   = fread(fid,1,'int16')';
dime.bitpix     = fread(fid,1,'int16')';
dime.dim_un0    = fread(fid,1,'int16')';
dime.pixdim     = fread(fid,8,'float32')';
dime.vox_offset = fread(fid,1,'float32')';
dime.roi_scale  = fread(fid,1,'float32')';
dime.funused1   = fread(fid,1,'float32')';
dime.funused2   = fread(fid,1,'float32')';
dime.cal_max    = fread(fid,1,'float32')';
dime.cal_min    = fread(fid,1,'float32')';
dime.compressed = fread(fid,1,'int32')';
dime.verified   = fread(fid,1,'int32')';
dime.glmax      = fread(fid,1,'int32')';
dime.glmin      = fread(fid,1,'int32')';

return					% image_dimension
end

%---------------------------------------------------------------------
function [ hist ] = data_history_untouch0_LOADUNTOUCHNIIGZ(fid)

%struct data_history
%       {                                /* off + size      */
%       char descrip[80];                /* 0 + 80          */
%       char aux_file[24];               /* 80 + 24         */
%       char orient;                     /* 104 + 1         */
%       char originator[10];             /* 105 + 10        */
%       char generated[10];              /* 115 + 10        */
%       char scannum[10];                /* 125 + 10        */
%       char patient_id[10];             /* 135 + 10        */
%       char exp_date[10];               /* 145 + 10        */
%       char exp_time[10];               /* 155 + 10        */
%       char hist_un0[3];                /* 165 + 3         */
%       int views                        /* 168 + 4         */
%       int vols_added;                  /* 172 + 4         */
%       int start_field;                 /* 176 + 4         */
%       int field_skip;                  /* 180 + 4         */
%       int omax;                        /* 184 + 4         */
%       int omin;                        /* 188 + 4         */
%       int smax;                        /* 192 + 4         */
%       int smin;                        /* 196 + 4         */
%       };                               /* total=200 bytes */

v6 = version;
if str2num(v6(1))<6
   directchar = '*char';
else
   directchar = 'uchar=>char';
end

hist.descrip     = deblank(fread(fid,80,directchar)');
hist.aux_file    = deblank(fread(fid,24,directchar)');
hist.orient      = fread(fid, 1,'char')';
hist.originator  = fread(fid, 5,'int16')';
hist.generated   = deblank(fread(fid,10,directchar)');
hist.scannum     = deblank(fread(fid,10,directchar)');
hist.patient_id  = deblank(fread(fid,10,directchar)');
hist.exp_date    = deblank(fread(fid,10,directchar)');
hist.exp_time    = deblank(fread(fid,10,directchar)');
hist.hist_un0    = deblank(fread(fid, 3,directchar)');
hist.views       = fread(fid, 1,'int32')';
hist.vols_added  = fread(fid, 1,'int32')';
hist.start_field = fread(fid, 1,'int32')';
hist.field_skip  = fread(fid, 1,'int32')';
hist.omax        = fread(fid, 1,'int32')';
hist.omin        = fread(fid, 1,'int32')';
hist.smax        = fread(fid, 1,'int32')';
hist.smin        = fread(fid, 1,'int32')';

return					% data_history
end


function hdr = load_untouch_nii_hdr_LOADUNTOUCHNIIGZ(fileprefix, machine, filetype)

if filetype == 2
   fn = sprintf('%s.nii',fileprefix);
   
   if ~exist(fn)
      msg = sprintf('Cannot find file "%s.nii".', fileprefix);
      error(msg);
   end
else
   fn = sprintf('%s.hdr',fileprefix);
   
   if ~exist(fn)
      msg = sprintf('Cannot find file "%s.hdr".', fileprefix);
      error(msg);
   end
end

fid = fopen(fn,'r',machine);

if fid < 0,
   msg = sprintf('Cannot open file %s.',fn);
   error(msg);
else
   fseek(fid,0,'bof');
   hdr = read_header_untouch_LOADUNTOUCHNIIGZ(fid);
   fclose(fid);
end

return					% load_nii_hdr
end

%---------------------------------------------------------------------
function [ dsr ] = read_header_untouch_LOADUNTOUCHNIIGZ(fid)

%  Original header structures
%  struct dsr
%       {
%       struct header_key hk;            /*   0 +  40       */
%       struct image_dimension dime;     /*  40 + 108       */
%       struct data_history hist;        /* 148 + 200       */
%       };                               /* total= 348 bytes*/

dsr.hk   = header_key_untouch_LOADUNTOUCHNIIGZ(fid);
dsr.dime = image_dimension_untouch_LOADUNTOUCHNIIGZ(fid);
dsr.hist = data_history_untouch_LOADUNTOUCHNIIGZ(fid);

%  For Analyze data format
%
if ~strcmp(dsr.hist.magic, 'n+1') & ~strcmp(dsr.hist.magic, 'ni1')
   dsr.hist.qform_code = 0;
   dsr.hist.sform_code = 0;
end

return					% read_header
end

%---------------------------------------------------------------------
function [ hk ] = header_key_untouch_LOADUNTOUCHNIIGZ(fid)

fseek(fid,0,'bof');

%  Original header structures
%  struct header_key                     /* header key      */
%       {                                /* off + size      */
%       int sizeof_hdr                   /*  0 +  4         */
%       char data_type[10];              /*  4 + 10         */
%       char db_name[18];                /* 14 + 18         */
%       int extents;                     /* 32 +  4         */
%       short int session_error;         /* 36 +  2         */
%       char regular;                    /* 38 +  1         */
%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
%       };                               /* total=40 bytes  */
%
% int sizeof_header   Should be 348.
% char regular        Must be 'r' to indicate that all images and
%                     volumes are the same size.

v6 = version;
if str2num(v6(1))<6
   directchar = '*char';
else
   directchar = 'uchar=>char';
end

hk.sizeof_hdr    = fread(fid, 1,'int32')';	% should be 348!
hk.data_type     = deblank(fread(fid,10,directchar)');
hk.db_name       = deblank(fread(fid,18,directchar)');
hk.extents       = fread(fid, 1,'int32')';
hk.session_error = fread(fid, 1,'int16')';
hk.regular       = fread(fid, 1,directchar)';
hk.dim_info      = fread(fid, 1,'uchar')';

return					% header_key
end

%---------------------------------------------------------------------
function [ dime ] = image_dimension_untouch_LOADUNTOUCHNIIGZ(fid)

%  Original header structures
%  struct image_dimension
%       {                                /* off + size      */
%       short int dim[8];                /* 0 + 16          */
%       /*
%           dim[0]      Number of dimensions in database; usually 4.
%           dim[1]      Image X dimension;  number of *pixels* in an image row.
%           dim[2]      Image Y dimension;  number of *pixel rows* in slice.
%           dim[3]      Volume Z dimension; number of *slices* in a volume.
%           dim[4]      Time points; number of volumes in database
%       */
%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
%       short int intent_code;   % short int unused1;   /* 28 + 2 */
%       short int datatype;              /* 30 + 2          */
%       short int bitpix;                /* 32 + 2          */
%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
%       float pixdim[8];                 /* 36 + 32         */
%	/*
%		pixdim[] specifies the voxel dimensions:
%		pixdim[1] - voxel width, mm
%		pixdim[2] - voxel height, mm
%		pixdim[3] - slice thickness, mm
%		pixdim[4] - volume timing, in msec
%					..etc
%	*/
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

dime.dim        = fread(fid,8,'int16')';
dime.intent_p1  = fread(fid,1,'float32')';
dime.intent_p2  = fread(fid,1,'float32')';
dime.intent_p3  = fread(fid,1,'float32')';
dime.intent_code = fread(fid,1,'int16')';
dime.datatype   = fread(fid,1,'int16')';
dime.bitpix     = fread(fid,1,'int16')';
dime.slice_start = fread(fid,1,'int16')';
dime.pixdim     = fread(fid,8,'float32')';
dime.vox_offset = fread(fid,1,'float32')';
dime.scl_slope  = fread(fid,1,'float32')';
dime.scl_inter  = fread(fid,1,'float32')';
dime.slice_end  = fread(fid,1,'int16')';
dime.slice_code = fread(fid,1,'uchar')';
dime.xyzt_units = fread(fid,1,'uchar')';
dime.cal_max    = fread(fid,1,'float32')';
dime.cal_min    = fread(fid,1,'float32')';
dime.slice_duration = fread(fid,1,'float32')';
dime.toffset    = fread(fid,1,'float32')';
dime.glmax      = fread(fid,1,'int32')';
dime.glmin      = fread(fid,1,'int32')';

return					% image_dimension
end

%---------------------------------------------------------------------
function [ hist ] = data_history_untouch_LOADUNTOUCHNIIGZ(fid)

%  Original header structures
%  struct data_history
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

v6 = version;
if str2num(v6(1))<6
   directchar = '*char';
else
   directchar = 'uchar=>char';
end

hist.descrip     = deblank(fread(fid,80,directchar)');
hist.aux_file    = deblank(fread(fid,24,directchar)');
hist.qform_code  = fread(fid,1,'int16')';
hist.sform_code  = fread(fid,1,'int16')';
hist.quatern_b   = fread(fid,1,'float32')';
hist.quatern_c   = fread(fid,1,'float32')';
hist.quatern_d   = fread(fid,1,'float32')';
hist.qoffset_x   = fread(fid,1,'float32')';
hist.qoffset_y   = fread(fid,1,'float32')';
hist.qoffset_z   = fread(fid,1,'float32')';
hist.srow_x      = fread(fid,4,'float32')';
hist.srow_y      = fread(fid,4,'float32')';
hist.srow_z      = fread(fid,4,'float32')';
hist.intent_name = deblank(fread(fid,16,directchar)');
hist.magic       = deblank(fread(fid,4,directchar)');

return					% data_history
end

function [img,hdr] = load_untouch_nii_img_LOADUNTOUCHNIIGZ(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB,slice_idx)

if ~exist('hdr','var') | ~exist('filetype','var') | ~exist('fileprefix','var') | ~exist('machine','var')
   error('Usage: [img,hdr] = load_nii_img(hdr,filetype,fileprefix,machine,[img_idx],[dim5_idx],[dim6_idx],[dim7_idx],[old_RGB],[slice_idx]);');
end

if ~exist('img_idx','var') | isempty(img_idx) | hdr.dime.dim(5)<1
   img_idx = [];
end

if ~exist('dim5_idx','var') | isempty(dim5_idx) | hdr.dime.dim(6)<1
   dim5_idx = [];
end

if ~exist('dim6_idx','var') | isempty(dim6_idx) | hdr.dime.dim(7)<1
   dim6_idx = [];
end

if ~exist('dim7_idx','var') | isempty(dim7_idx) | hdr.dime.dim(8)<1
   dim7_idx = [];
end

if ~exist('old_RGB','var') | isempty(old_RGB)
   old_RGB = 0;
end

if ~exist('slice_idx','var') | isempty(slice_idx) | hdr.dime.dim(4)<1
   slice_idx = [];
end

%  check img_idx
%
if ~isempty(img_idx) & ~isnumeric(img_idx)
   error('"img_idx" should be a numerical array.');
end

if length(unique(img_idx)) ~= length(img_idx)
   error('Duplicate image index in "img_idx"');
end

if ~isempty(img_idx) & (min(img_idx) < 1 | max(img_idx) > hdr.dime.dim(5))
   max_range = hdr.dime.dim(5);
   
   if max_range == 1
      error(['"img_idx" should be 1.']);
   else
      range = ['1 ' num2str(max_range)];
      error(['"img_idx" should be an integer within the range of [' range '].']);
   end
end

%  check dim5_idx
%
if ~isempty(dim5_idx) & ~isnumeric(dim5_idx)
   error('"dim5_idx" should be a numerical array.');
end

if length(unique(dim5_idx)) ~= length(dim5_idx)
   error('Duplicate index in "dim5_idx"');
end

if ~isempty(dim5_idx) & (min(dim5_idx) < 1 | max(dim5_idx) > hdr.dime.dim(6))
   max_range = hdr.dime.dim(6);
   
   if max_range == 1
      error(['"dim5_idx" should be 1.']);
   else
      range = ['1 ' num2str(max_range)];
      error(['"dim5_idx" should be an integer within the range of [' range '].']);
   end
end

%  check dim6_idx
%
if ~isempty(dim6_idx) & ~isnumeric(dim6_idx)
   error('"dim6_idx" should be a numerical array.');
end

if length(unique(dim6_idx)) ~= length(dim6_idx)
   error('Duplicate index in "dim6_idx"');
end

if ~isempty(dim6_idx) & (min(dim6_idx) < 1 | max(dim6_idx) > hdr.dime.dim(7))
   max_range = hdr.dime.dim(7);
   
   if max_range == 1
      error(['"dim6_idx" should be 1.']);
   else
      range = ['1 ' num2str(max_range)];
      error(['"dim6_idx" should be an integer within the range of [' range '].']);
   end
end

%  check dim7_idx
%
if ~isempty(dim7_idx) & ~isnumeric(dim7_idx)
   error('"dim7_idx" should be a numerical array.');
end

if length(unique(dim7_idx)) ~= length(dim7_idx)
   error('Duplicate index in "dim7_idx"');
end

if ~isempty(dim7_idx) & (min(dim7_idx) < 1 | max(dim7_idx) > hdr.dime.dim(8))
   max_range = hdr.dime.dim(8);
   
   if max_range == 1
      error(['"dim7_idx" should be 1.']);
   else
      range = ['1 ' num2str(max_range)];
      error(['"dim7_idx" should be an integer within the range of [' range '].']);
   end
end

%  check slice_idx
%
if ~isempty(slice_idx) & ~isnumeric(slice_idx)
   error('"slice_idx" should be a numerical array.');
end

if length(unique(slice_idx)) ~= length(slice_idx)
   error('Duplicate index in "slice_idx"');
end

if ~isempty(slice_idx) & (min(slice_idx) < 1 | max(slice_idx) > hdr.dime.dim(4))
   max_range = hdr.dime.dim(4);
   
   if max_range == 1
      error(['"slice_idx" should be 1.']);
   else
      range = ['1 ' num2str(max_range)];
      error(['"slice_idx" should be an integer within the range of [' range '].']);
   end
end

[img,hdr] = read_image_LOADUNTOUCHNIIGZ(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB,slice_idx);

return					% load_nii_img
end

%---------------------------------------------------------------------
function [img,hdr] = read_image_LOADUNTOUCHNIIGZ(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB,slice_idx)

switch filetype
   case {0, 1}
      fn = [fileprefix '.img'];
   case 2
      fn = [fileprefix '.nii'];
end

fid = fopen(fn,'r',machine);

if fid < 0,
   msg = sprintf('Cannot open file %s.',fn);
   error(msg);
end

%  Set bitpix according to datatype
%
%  /*Acceptable values for datatype are*/
%
%     0 None                     (Unknown bit per voxel) % DT_NONE, DT_UNKNOWN
%     1 Binary                         (ubit1, bitpix=1) % DT_BINARY
%     2 Unsigned char         (uchar or uint8, bitpix=8) % DT_UINT8, NIFTI_TYPE_UINT8
%     4 Signed short                  (int16, bitpix=16) % DT_INT16, NIFTI_TYPE_INT16
%     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32
%    16 Floating point    (single or float32, bitpix=32) % DT_FLOAT32, NIFTI_TYPE_FLOAT32
%    32 Complex, 2 float32      (Use float32, bitpix=64) % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
%    64 Double precision  (double or float64, bitpix=64) % DT_FLOAT64, NIFTI_TYPE_FLOAT64
%   128 uint8 RGB                 (Use uint8, bitpix=24) % DT_RGB24, NIFTI_TYPE_RGB24
%   256 Signed char            (schar or int8, bitpix=8) % DT_INT8, NIFTI_TYPE_INT8
%   511 Single RGB              (Use float32, bitpix=96) % DT_RGB96, NIFTI_TYPE_RGB96
%   512 Unsigned short               (uint16, bitpix=16) % DT_UNINT16, NIFTI_TYPE_UNINT16
%   768 Unsigned integer             (uint32, bitpix=32) % DT_UNINT32, NIFTI_TYPE_UNINT32
%  1024 Signed long long              (int64, bitpix=64) % DT_INT64, NIFTI_TYPE_INT64
%  1280 Unsigned long long           (uint64, bitpix=64) % DT_UINT64, NIFTI_TYPE_UINT64
%  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128
%  1792 Complex128, 2 float64  (Use float64, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128
%  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128
%
switch hdr.dime.datatype
   case   1,
      hdr.dime.bitpix = 1;  precision = 'ubit1';
   case   2,
      hdr.dime.bitpix = 8;  precision = 'uint8';
   case   4,
      hdr.dime.bitpix = 16; precision = 'int16';
   case   8,
      hdr.dime.bitpix = 32; precision = 'int32';
   case  16,
      hdr.dime.bitpix = 32; precision = 'float32';
   case  32,
      hdr.dime.bitpix = 64; precision = 'float32';
   case  64,
      hdr.dime.bitpix = 64; precision = 'float64';
   case 128,
      hdr.dime.bitpix = 24; precision = 'uint8';
   case 256
      hdr.dime.bitpix = 8;  precision = 'int8';
   case 511
      hdr.dime.bitpix = 96; precision = 'float32';
   case 512
      hdr.dime.bitpix = 16; precision = 'uint16';
   case 768
      hdr.dime.bitpix = 32; precision = 'uint32';
   case 1024
      hdr.dime.bitpix = 64; precision = 'int64';
   case 1280
      hdr.dime.bitpix = 64; precision = 'uint64';
   case 1792,
      hdr.dime.bitpix = 128; precision = 'float64';
   otherwise
      error('This datatype is not supported');
end

tmp = hdr.dime.dim(2:end);
tmp(find(tmp < 1)) = 1;
hdr.dime.dim(2:end) = tmp;

%  move pointer to the start of image block
%
switch filetype
   case {0, 1}
      fseek(fid, 0, 'bof');
   case 2
      fseek(fid, hdr.dime.vox_offset, 'bof');
end

%  Load whole image block for old Analyze format or binary image;
%  otherwise, load images that are specified in img_idx, dim5_idx,
%  dim6_idx, and dim7_idx
%
%  For binary image, we have to read all because pos can not be
%  seeked in bit and can not be calculated the way below.
%
if hdr.dime.datatype == 1 | isequal(hdr.dime.dim(4:8),ones(1,5)) | ...
      (isempty(img_idx) & isempty(dim5_idx) & isempty(dim6_idx) & isempty(dim7_idx) & isempty(slice_idx))
   
   %  For each frame, precision of value will be read
   %  in img_siz times, where img_siz is only the
   %  dimension size of an image, not the byte storage
   %  size of an image.
   %
   img_siz = prod(hdr.dime.dim(2:8));
   
   %  For complex float32 or complex float64, voxel values
   %  include [real, imag]
   %
   if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
      img_siz = img_siz * 2;
   end
   
   %MPH: For RGB24, voxel values include 3 separate color planes
   %
   if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
      img_siz = img_siz * 3;
   end
   
   img = fread(fid, img_siz, sprintf('*%s',precision));
   
   d1 = hdr.dime.dim(2);
   d2 = hdr.dime.dim(3);
   d3 = hdr.dime.dim(4);
   d4 = hdr.dime.dim(5);
   d5 = hdr.dime.dim(6);
   d6 = hdr.dime.dim(7);
   d7 = hdr.dime.dim(8);
   
   if isempty(slice_idx)
      slice_idx = 1:d3;
   end
   
   if isempty(img_idx)
      img_idx = 1:d4;
   end
   
   if isempty(dim5_idx)
      dim5_idx = 1:d5;
   end
   
   if isempty(dim6_idx)
      dim6_idx = 1:d6;
   end
   
   if isempty(dim7_idx)
      dim7_idx = 1:d7;
   end
else
   
   d1 = hdr.dime.dim(2);
   d2 = hdr.dime.dim(3);
   d3 = hdr.dime.dim(4);
   d4 = hdr.dime.dim(5);
   d5 = hdr.dime.dim(6);
   d6 = hdr.dime.dim(7);
   d7 = hdr.dime.dim(8);
   
   if isempty(slice_idx)
      slice_idx = 1:d3;
   end
   
   if isempty(img_idx)
      img_idx = 1:d4;
   end
   
   if isempty(dim5_idx)
      dim5_idx = 1:d5;
   end
   
   if isempty(dim6_idx)
      dim6_idx = 1:d6;
   end
   
   if isempty(dim7_idx)
      dim7_idx = 1:d7;
   end
   
   %ROMAN: begin
   roman = 1;
   if(roman)
      
      %  compute size of one slice
      %
      img_siz = prod(hdr.dime.dim(2:3));
      
      %  For complex float32 or complex float64, voxel values
      %  include [real, imag]
      %
      if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
         img_siz = img_siz * 2;
      end
      
      %MPH: For RGB24, voxel values include 3 separate color planes
      %
      if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
         img_siz = img_siz * 3;
      end
      
      % preallocate img
      img = zeros(img_siz, length(slice_idx)*length(img_idx)*length(dim5_idx)*length(dim6_idx)*length(dim7_idx) );
      currentIndex = 1;
   else
      img = [];
   end; %if(roman)
   % ROMAN: end
   
   for i7=1:length(dim7_idx)
      for i6=1:length(dim6_idx)
         for i5=1:length(dim5_idx)
            for t=1:length(img_idx)
               for s=1:length(slice_idx)
                  
                  %  Position is seeked in bytes. To convert dimension size
                  %  to byte storage size, hdr.dime.bitpix/8 will be
                  %  applied.
                  %
                  pos = sub2ind([d1 d2 d3 d4 d5 d6 d7], 1, 1, slice_idx(s), ...
                     img_idx(t), dim5_idx(i5),dim6_idx(i6),dim7_idx(i7)) -1;
                  pos = pos * hdr.dime.bitpix/8;
                  
                  % ROMAN: begin
                  if(roman)
                     % do nothing
                  else
                     img_siz = prod(hdr.dime.dim(2:3));
                     
                     %  For complex float32 or complex float64, voxel values
                     %  include [real, imag]
                     %
                     if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
                        img_siz = img_siz * 2;
                     end
                     
                     %MPH: For RGB24, voxel values include 3 separate color planes
                     %
                     if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
                        img_siz = img_siz * 3;
                     end
                  end; % if (roman)
                  % ROMAN: end
                  
                  if filetype == 2
                     fseek(fid, pos + hdr.dime.vox_offset, 'bof');
                  else
                     fseek(fid, pos, 'bof');
                  end
                  
                  %  For each frame, fread will read precision of value
                  %  in img_siz times
                  %
                  % ROMAN: begin
                  if(roman)
                     img(:,currentIndex) = fread(fid, img_siz, sprintf('*%s',precision));
                     currentIndex = currentIndex +1;
                  else
                     img = [img fread(fid, img_siz, sprintf('*%s',precision))];
                  end; %if(roman)
                  % ROMAN: end
                  
               end
            end
         end
      end
   end
end

%  For complex float32 or complex float64, voxel values
%  include [real, imag]
%
if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
   img = reshape(img, [2, length(img)/2]);
   img = complex(img(1,:)', img(2,:)');
end

fclose(fid);

%  Update the global min and max values
%
hdr.dime.glmax = double(max(img(:)));
hdr.dime.glmin = double(min(img(:)));

%  old_RGB treat RGB slice by slice, now it is treated voxel by voxel
%
if old_RGB & hdr.dime.datatype == 128 & hdr.dime.bitpix == 24
   % remove squeeze
   img = (reshape(img, [hdr.dime.dim(2:3) 3 length(slice_idx) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
   img = permute(img, [1 2 4 3 5 6 7 8]);
elseif hdr.dime.datatype == 128 & hdr.dime.bitpix == 24
   % remove squeeze
   img = (reshape(img, [3 hdr.dime.dim(2:3) length(slice_idx) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
   img = permute(img, [2 3 4 1 5 6 7 8]);
elseif hdr.dime.datatype == 511 & hdr.dime.bitpix == 96
   img = double(img(:));
   img = single((img - min(img))/(max(img) - min(img)));
   % remove squeeze
   img = (reshape(img, [3 hdr.dime.dim(2:3) length(slice_idx) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
   img = permute(img, [2 3 4 1 5 6 7 8]);
else
   % remove squeeze
   img = (reshape(img, [hdr.dime.dim(2:3) length(slice_idx) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
end

if ~isempty(slice_idx)
   hdr.dime.dim(4) = length(slice_idx);
end

if ~isempty(img_idx)
   hdr.dime.dim(5) = length(img_idx);
end

if ~isempty(dim5_idx)
   hdr.dime.dim(6) = length(dim5_idx);
end

if ~isempty(dim6_idx)
   hdr.dime.dim(7) = length(dim6_idx);
end

if ~isempty(dim7_idx)
   hdr.dime.dim(8) = length(dim7_idx);
end

return						% read_image
end

