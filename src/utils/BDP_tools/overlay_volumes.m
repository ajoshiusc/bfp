% 
% BDP BrainSuite Diffusion Pipeline
% 
% Copyright (C) 2019 The Regents of the University of California and
% the University of Southern California
% 
% Created by Chitresh Bhushan, Divya Varadarajan, Justin P. Haldar, Anand A. Joshi,
%            David W. Shattuck, and Richard M. Leahy
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
% 


function Img = overlay_volumes(vol1, vol2, varargin)
% VOL1 and VOL2 should be 3D matrix & the intensity should be
% normalized in range [0 1] or set CLIM to normalize intensity. 
%
% CLIM (optional) can be used to normalize the intensity. When used, it must NOT come 
%      immediately after TYPE argument. 
%
% TYPE (optional) can be 'rview', ['redgreen'], 'greenblue', 'yellowblue'
%     'rview' type overlay may hide low intensity areas. 
%
% TYPE_CLIM (optional) is percentile range for min-max intensity. When used, it must 
%      come immediately after TYPE argument. Default value of TYPE_CLIM is [70 94]. 
%      See normalize_intensity.m 
% 
% Usage: 
%       overlay_volumes(vol1, vol2)
%       overlay_volumes(vol1, vol2, dim)
%       overlay_volumes(vol1, vol2, clim)
%       overlay_volumes(vol1, vol2, type)
%       overlay_volumes(vol1, vol2, type, type_clim)
%       overlay_volumes(vol1, vol2, ...)
%

vol1 = double(vol1);
vol2 = double(vol2);

for n = 3:nargin
   k = n-2;
   if ischar(varargin{k})
      type = varargin{k};
   elseif isnumeric(varargin{k}) && length(varargin{k}) == 2
      if k>1 && ischar(varargin{k-1})
         type_clim = varargin{k};
      else
         clim = varargin{k};
      end
   else
      dim = varargin{k};
   end
end

if ~exist('clim', 'var')
   clim = [0 1];
end

if ~exist('type', 'var')
   type = 'redgreen';
end

if ~exist('type_clim', 'var')   
   type_clim = [70 94];
end

if ~exist('dim', 'var')
   dim = 3;
end

vol1 = (vol1-clim(1))/(clim(2)-clim(1));
vol1(vol1<0) = 0;
vol1(vol1>1) = 1;

vol2 = (vol2-clim(1))/(clim(2)-clim(1));
vol2(vol2<0) = 0;
vol2(vol2>1) = 1;

Img = zeros([size(vol1) 3]);

switch type
   case 'redgreen'
      Img(:,:,:,1) = vol1;
      Img(:,:,:,2) = vol2;
      Img(:,:,:,3) = 0;
      
   case 'rview'
      alpha1 = 0.5;
      alpha2 = 1;
      alpha_n = alpha1 + alpha2*(1-alpha1);
      p = normalize_intensity(vol2, type_clim);
      
      Img(:,:,:,1) = vol1*alpha1;
      Img(:,:,:,3) = vol1*alpha1;
      Img(:,:,:,2) = vol1*alpha1 + p*alpha2*(1-alpha1);
      Img = Img/alpha_n;
      
   case 'edge'
      %       alpha1 = 0.7;
      %       alpha2 = 1;
      %       alpha_n = alpha1 + alpha2*(1-alpha1);
      %
      %       Img(:,:,:,2) = vol1*alpha1;
      %       Img(:,:,:,1) = vol1*alpha1;
      %       Img(:,:,:,3) = vol1*alpha1 + vol2*alpha2*(1-alpha1);
      %       Img = Img/alpha_n;
      
      vol1(vol2>0) = 0;
      Img(:,:,:,2) = vol1;
      Img(:,:,:,3) = vol1;
      Img(:,:,:,1) = vol1 + vol2;

   case 'greenblue'
      Img(:,:,:,2) = vol1;
      Img(:,:,:,3) = vol2;
      Img(:,:,:,1) = 0;
      
   case 'yellowblue'
      Img(:,:,:,1) = vol1;
      Img(:,:,:,2) = vol1;
      Img(:,:,:,3) = vol2;
      
   otherwise
      error('unknown type')
end
clear vol1 vol2


if nargout<1   
   switch dim
      case 3
         E = permute(Img, [2 1 4 3]);
      case 2
         E = permute(Img, [3 1 4 2]);
      case 1
         E = permute(Img, [3 2 4 1]);
      otherwise
         error('dim can not be more than 3');
   end
   E = flipdim(E, 1);
   figure;
   
   imdisp(E, 'Size', 1);
   clear Img
end


end
