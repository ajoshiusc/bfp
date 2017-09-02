% SVReg: Surface-Constrained Volumetric Registration
% Copyright (C) 2017 The Regents of the University of California and the University of Southern California
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


function [nii, reorient_matrix, sform_new, niifileFlag] = load_nii_BIG_Lab(fname)
% Wrapper around many custom functions to ease life. It does following
% while loading nifti-1 file : 
%   1. Adds sform matrix, if sform does not exist.
%   2. Applies reorientation to volume to make sform cannonical (sform & other nifti headers are
%      also appropriately corrected).
%   3. Updates qform to reflect new reoriented sform.
%
% Details: 
%   fname - filename of nifti file; .nii or .nii.gz
%   nii - loaded nifti structure
%   reorient_matrix, sform_new, niifileFlag - see corresponding functions for details. 
%
% The loaded nifti structure can be saved to disk by calling either of following: 
%   save_untouch_nii_gz(nii, 'output_filename.nii.gz')
%   save_untouch_nii(nii, 'output_filename.nii') % using Jimmy Shen's toolbox
%
% Note that this function: 
%   1. Prefers sform over qform, when both are present (can be modified easily, if required)
%   2. Uses some intermediate file I/O for some operation - so could be a little slower
%

workdir = tempname(); mkdir(workdir);
temp_fname = fullfile(workdir, [Random_String(16) '.nii.gz']);

[niifileFlag, outFile] = check_nifti_file(fname, workdir);
[nii_out, reorient_matrix, sform_new] = reorient_nifti_sform(outFile, temp_fname);
nii = add_qform(nii_out);

rmdir(workdir, 's');

end
