function M = par2affineMat(t,r,s,h)
% Creates a 3D affine transformation matrix from translation (t),
% rotation(r), scale(s) and shear(h) parameters.
%
%   t: [t_x t_y t_z]
%   r: [rotateX rotateY rotateZ] in degrees
%   s: [scaleX scaleY scaleZ]
%   h: [ShearXY, ShearXZ, ShearYX, ShearYZ, ShearZX, ShearZY]
%

if ~exist('r','var') || isempty(r)
   r = [0 0 0];
end

if ~exist('s','var') || isempty(s)
   s = [1 1 1];
end

if ~exist('h','var') || isempty(h)
   h = [0 0 0 0 0 0];
end

M = mat_tra(t) * mat_scale(s) * mat_rot(r) * mat_shear(h);

end

function M = mat_rot(r)
r = r * (pi/180);
Rx = ...
   [1 0 0 0;
   0 cos(r(1)) -sin(r(1)) 0;
   0 sin(r(1)) cos(r(1)) 0;
   0 0 0 1];

Ry = ...
   [cos(r(2)) 0 sin(r(2)) 0;
   0 1 0 0;
   -sin(r(2)) 0 cos(r(2)) 0;
   0 0 0 1];

Rz = ...
   [cos(r(3)) -sin(r(3)) 0 0;
   sin(r(3)) cos(r(3)) 0 0;
   0 0 1 0;
   0 0 0 1];

M = Rx*Ry*Rz;
end

function M = mat_scale(s)
M = ...
   [s(1) 0    0    0;
   0    s(2) 0    0;
   0    0    s(3) 0;
   0    0    0    1];
end

function M = mat_shear(h)
M = ...
   [1    h(1) h(2) 0;
   h(3) 1    h(4) 0;
   h(5) h(6) 1    0;
   0 0 0 1];
end

function M = mat_tra(t)
M = ...
   [1 0 0 t(1);
   0 1 0 t(2);
   0 0 1 t(3);
   0 0 0 1];
end
