function [A, x_out, y_out, z_out] = my_imresize3d(V, scale, tsize, ntype)
% This function resizes a 3D image volume to new dimensions. More accurate computation of image grid
% for output image. 
%   V: The input image volume
%   scale: scalar scaling factor - could be a scalar or vector of size 3; when used set tsize to [];
%   tsize: new dimensions -vector of size 3; when used set scale to [];
%   ntype: (optional) Type of interpolation ('nearest', 'linear', or 'cubic') - default: linear
%
% Usage: 
%   A = my_imresize3d(V, [], [20 30 45])
%   A = my_imresize3d(V, 0.85, [])
%   A = my_imresize3d(V, [0.5 0.6 2], [])
%   A = my_imresize3d(..., 'cubic')
% 
% Testing: 
%    sz = size(sImgP);
%    tsz = round(size(sImgP)*2.5);
%    temp1 = my_imresize3d(sImgP, [], tsz, 'linear'); % upsample
%    temp2 = my_imresize3d(temp1, [], sz, 'linear');  % downsample
%    overlay_volumes(sImgP, temp2)
%    display_volume(temp2-sImgP)
%    
%    % test 2
%    testimg = my_imresize3d(sImgP, [], 4*sz, 'linear');
%    tsz = round(4*size(sImgP)/2.85);
%    temp1 = my_imresize3d(testimg, [], tsz, 'linear'); % downsample
%    temp2 = my_imresize3d(temp1, [], 4*sz, 'linear');  % upsample
%    overlay_volumes(testimg, temp2)
%

% Check the inputs
if ~isempty(scale) && ~isempty(tsize)
   error('Only one of scale & tsize should be set. Other one must be empty.')
end

if ~exist('ntype', 'var')
   ntype = 'linear';
end

vsize = size(V);
[x_out,y_out,z_out] = ndgrid(1:vsize(1), 1:vsize(2), 1:vsize(3)); % in case no interpolation is required
if ~isempty(scale)
   if scale==1
      A = V;
      return;
   else
      tsize = round(vsize.*scale);
   end
end

if ~isempty(tsize) && isequal(tsize, vsize)
   A = V;
   return;
end

delta = vsize./tsize;
offset = delta/2 - 0.5; % -0.5 to start indexing from 0 (& 0.5 = voxel extend from voxel-center wrt V coord)
                        % delta/2 to shift to center of voxel in target coord
ulim = (tsize-1).*delta;
[x,y,z] = ndgrid(0:(vsize(1)-1), 0:(vsize(2)-1), 0:(vsize(3)-1));
[x_out,y_out,z_out] = ndgrid(0:delta(1):ulim(1), 0:delta(2):ulim(2), 0:delta(3):ulim(3));
x_out = x_out + offset(1);
y_out = y_out + offset(2);
z_out = z_out + offset(3);

A = interpnNNboundary(x,y,z, double(V), x_out, y_out, z_out, ntype, 0, 0.15);

% convert to indexing from 1 (matlab style)
x_out = x_out + 1;
y_out = y_out + 1;
z_out = z_out + 1;

end

