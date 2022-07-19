% SVReg: Surface-Constrained Volumetric Registration
% Copyright (C) 2019 The Regents of the University of California and the University of Southern California
% Created by Anand A. Joshi, Chitresh Bhushan, David W. Shattuck, Richard M. Leahy 
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; version 2.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
% USA.

function svreg_resample(varargin)
%
% Description:
% This function resamples 3D or 4D input NIFTI-1 file.
%
% Usage:
% svreg_resample(in_file,out_file,[dim],[dx],[dy],[dz],[method],[extrapval])
% while first two arguments are required, others are [OPTIONAL]
%
% Arguments:
% in_file: input nifti file in nii.gz file format
% out_file: output nifti file in nii.gz format
% [dim]: if this is '-size' output size is specified
%        if this is '-res' output pixel resolution is specified
%       default is -res 1 1 1 (istotropic sampling)
% [dx],[dy],[dz]: if the dim='-size', then this is interpreted as output size
%                 if the dim=-res', then this is interpred as voxel resolution in x,y,z directions
% [method]: interpolation methods: 'linear', 'nearest','cubic','spline' the
% default is 'linear'.
% [extrapval]: Values for extrapolation outside the grid on the boundary.
% the default is median of all corners.
%
% If the functions is called with only input and output file names, it
% resamples to isotropic (1mm cubic) resolution

p = inputParser;
defaultDx = '1';   defaultDy = '1';   defaultDz = '1';
defaultDim = '-res';
defaultMethod='linear';
defaultExtapVal='';
addRequired(p,'infile',@isstr);
addRequired(p,'outfile',@isstr);
addOptional(p,'dim',defaultDim,@isstr);
addOptional(p,'dx',defaultDx,@(x) ischar(x)||isnumeric(x));
addOptional(p,'dy',defaultDy,@(x) ischar(x)||isnumeric(x));
addOptional(p,'dz',defaultDz,@(x) ischar(x)||isnumeric(x));
addOptional(p,'method',defaultMethod,@isstr);
addOptional(p,'extrapval',defaultExtapVal,@(x) ischar(x)||isnumeric(x));

% parse the input arguments
parse(p,varargin{:})
if ischar(p.Results.dx)
    dx=str2double(p.Results.dx);
    dy=str2double(p.Results.dy);
    dz=str2double(p.Results.dz);
elseif isnumeric(p.Results.dx)
    dx=p.Results.dx;
    dy=p.Results.dy;
    dz=p.Results.dz;
end
dim=p.Results.dim;
infile=p.Results.infile;outfile=p.Results.outfile;
mthd=p.Results.method;
extrapval=str2double(p.Results.extrapval);
v=load_nii_BIG_Lab(infile);

% parse dx dy dz
if strcmp(dim,'-size')
    SZ=[dx,dy,dz]; RES=[];
elseif strcmp(dim,'-res')
    SZ=[]; RES=[dx,dy,dz];
else
    error('dim parameter is wrong.');
end

%extrap value is median of values at all corners in the image

if isempty(extrapval)
    extrapval=median([v.img(1,1,1),v.img(end,1,1),v.img(1,end,1),...
        v.img(1,1,end),v.img(end,1,end),v.img(1,end,end),...
        v.img(end,end,1),v.img(end,end,end)]);
end
% perform the resampling
vr=resample_vol_res(v,RES,SZ,mthd,extrapval);


% Make the header look like BrainSuite header
vr.hdr.hist.srow_x(1:3)=[vr.hdr.dime.pixdim(2),0,0];
vr.hdr.hist.srow_y(1:3)=[0,vr.hdr.dime.pixdim(3),0];
vr.hdr.hist.srow_z(1:3)=[0,0,vr.hdr.dime.pixdim(4)];
vr.hdr.dime.scl_slope= 0;

% save the output

save_untouch_nii(vr,outfile);
