function VI = interpnNNboundary(X1,X2,X3, V, Y1,Y2,Y3, method, extrapval, nbr_dist_thr)
% NOT a generic function. Should be used carefully!
%    X1,X2,X3, MUST be ndgrid type grid starting from 0 to (size(V)-1). X1<0 and X1>(size(V)-1) are
%    outside the volume V. Also, only intended to work with 3D images (no 2D or >3D). 
%
% Interpolation wrapper around Matlab's interpn to interpolate voxels near boundaries which lie
% outside the volume because of numerical error etc. in X1, X2, X3. nbr_dist_thr defines the
% boundary threshold.
% 
% Usage: 
%   VI = interpnNNboundary(X1,X2,X3, V, Y1,Y2,Y3)
%   VI = interpnNNboundary(X1,X2,X3, V, Y1,Y2,Y3, method)
%   VI = interpnNNboundary(X1,X2,X3, V, Y1,Y2,Y3, method, extrapval)
%   VI = interpnNNboundary(X1,X2,X3, V, Y1,Y2,Y3, method, extrapval, nbr_dist_thr)

if nargin<7
   error('Not enough input. See usage')
   
elseif nargin==7
   method = 'linear';
   extrapval = 0;
   nbr_dist_thr = 0.1;
   
elseif nargin==8
   extrapval = 0;
   nbr_dist_thr = 0.1;
   
elseif nargin==9
   nbr_dist_thr = 0.1;
end

VI = interpn(X1, X2, X3, double(V), Y1, Y2, Y3, method);
nan_ind = find(isnan(VI));
VI(nan_ind) = extrapval;

if ~isempty(nan_ind)
   [bdr_ind, X_bdr, Y_bdr, Z_bdr] = nearest_interp_boundary(Y1, Y2, Y3, size(V), nan_ind, nbr_dist_thr);
   VI(bdr_ind) = interpn(X1, X2, X3, double(V), X_bdr, Y_bdr, Z_bdr, method); % 'method' b/c 'nearest' results looks wierd (although correct)
end

end


function [boundary_ind, x, y, z] = nearest_interp_boundary(Y1, Y2, Y3, vol_size, nan_ind, nbr_dist_thr)
% Inputs: grid locations. Original grid locations are (0:vol_size-1)
%         index of nan points (corresponding to matix of grid locations)
%         neighbor distance (city-block) threshold (which would be replaced by neighbors)
%
% Outputs: index of boundary points (corresponding to matix of grid locations)
%          locations of neighbors

x = Y1(nan_ind);
y = Y2(nan_ind);
z = Y3(nan_ind);
boundary_ind = nan_ind;

% throw away voxels far away from boundary
m = x>=-nbr_dist_thr & x<=(vol_size(1)-1 + nbr_dist_thr) ...
   & y>=-nbr_dist_thr & y<=(vol_size(2)-1 + nbr_dist_thr) ...
   & z>=-nbr_dist_thr & z<=(vol_size(3)-1 + nbr_dist_thr);
boundary_ind(~m) = [];
x(~m) = [];
y(~m) = [];
z(~m) = [];

x = min(max(x,0),vol_size(1)-1);
y = min(max(y,0),vol_size(2)-1);
z = min(max(z,0),vol_size(3)-1);
end
