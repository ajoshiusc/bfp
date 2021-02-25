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


function [cp_mask_indx, cp_col_index, cp_col_sten] = bspline_create_mask_index(Ogrid_size, Spacing, I_size)
% Creates index (of image voxels) to be included inside area of influence of each cubic b-spline
% control point. Note that maximum voxel-index could be (prod(I_size)+1), and in that case
% corresponding voxel-intensity should be treated as 0. This is intended to be used for faster
% computation of analytical gradients. See Method1 and Method2 below. 
%
% Voxel-indices for (i,j,k)th control-points can be extracted by:
%      ind = cp_mask_indx(:, cp_col_index(i,j,k)); 
% 
%
% Method 1 - Easier to understand 
% 
%    % d_P is the image on which masking is to applied
%    d_Prod = zeros(numel(d_P)+1, 1); % for fake 0-voxel at (prod(I_size)+1)th location
%    O_grad1 = zeros(Ogrid_size(1:3)); 
%    for zi=1:4
%       for zj=1:4
%          for zk=1:4
%
%             d_Prod(1:end-1) = d_P(:); % copy image -- can this be done any faster?
%             
%             [xx, yy, zz] = ndgrid(zi:4:Ogrid_size(1), zj:4:Ogrid_size(2), zk:4:Ogrid_size(3));
%             CP_ind = sub2ind(Ogrid_size(1:3), xx(:), yy(:), zz(:)); % control-points in the set
%             
%             d_msk = sum(d_Prod(cp_mask_indx(:, cp_col_index(CP_ind))), 1);             
%             O_grad1(CP_ind) = d_msk(:);
%          end
%       end
%    end
% 
% 
% Method 2 - Slightly faster than method 1
% 
%    % d_P is the image on which masking is to applied
%    d_Prod = zeros(numel(d_P)+1, 1); % for fake 0-voxel at (prod(I_size)+1)th location
%    O_grad2 = zeros(Ogrid_size(1:3)); 
%    for zi=1:4
%       for zj=1:4
%          for zk=1:4
%
%             d_Prod(1:end-1) = d_P(:); % copy image -- can this be done any faster?
%             
%             [xx, yy, zz] = ndgrid(zi:4:Ogrid_size(1), zj:4:Ogrid_size(2), zk:4:Ogrid_size(3));
%             CP_ind = sub2ind(Ogrid_size(1:3), xx(:), yy(:), zz(:)); % control-points in the set
%             
%             d_msk = sum(d_Prod(cp_mask_indx(:, cp_col_sten(zi,zj,zk,1):cp_col_sten(zi,zj,zk,2))), 1);             
%             O_grad2(CP_ind) = d_msk(:);
%          end
%       end
%    end
% 
% 



nVox = prod(I_size);
ind_sz = prod(Spacing(1:3)*4 + 1); % (max) number of voxels for each unique mask
cp_mask_indx = zeros(ind_sz, prod(Ogrid_size(1:3)), 'uint32'); % Each column contains voxel-indices
                                                         % corresponding to area of influence of
                                                         % the control-point

cp_col_index = zeros(Ogrid_size(1:3)); % for each control-point, stores the corresponding 
                                       % column number in cp_mask_indx

cp_col_sten = zeros(4,4,4,2); % for each set of control-points, stores the corresponding
                              % starting and ending column number in cp_mask_indx
                                       
cpc = 0; % control-point counter
for zk=1:4
   for zj=1:4
      for zi=1:4
         
         cp_col_sten(zi,zj,zk,1) = cpc + 1; % starting location
         
         for k = zk:4:Ogrid_size(3)
            for j = zj:4:Ogrid_size(2)
               for i = zi:4:Ogrid_size(1)
                  
                  cpc = cpc + 1;
                  cp_col_index(i, j, k) = cpc;
                  
                  ind = zeros(ind_sz,1) + (nVox+1); % defaults to index of fake voxel with zero value
                  loc = ([i j k]-1).*Spacing + 1; % location in image grid
                  
                  xmin = max(1, loc(1)-2*Spacing(1));
                  xmax = min(I_size(1), loc(1)+2*Spacing(1));
                  
                  ymin = max(1, loc(2)-2*Spacing(2));
                  ymax = min(I_size(2), loc(2)+2*Spacing(2));
                  
                  zmin = max(1, loc(3)-2*Spacing(3));
                  zmax = min(I_size(3), loc(3)+2*Spacing(3));
                  
                  [xx,yy,zz] = ndgrid(xmin:xmax, ymin:ymax, zmin:zmax);
                  ind(1:numel(xx)) = sub2ind(I_size, xx(:), yy(:), zz(:));
                  
                  cp_mask_indx(:, cpc) = ind;
               end
            end
         end
         cp_col_sten(zi,zj,zk,2) = cpc; % ending location
      end
   end
end
end
