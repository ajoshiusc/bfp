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


function x_smooth = smoothGaussian1D(x, sgma)
% Replacement for smooth() and normcdf() to get rid of some toolbox dependencies. 
%
% Smooths a vector x using convolution with Gaussian pdf. End points are repeated for padding. 
%

xt = 0:6*sgma;
xt = [-1*xt(end:-1:2) xt]; % always odd length

if length(xt)<2
   x_smooth = x;
   
else
   P = normcdfBDP(xt, 0, sgma);
   u = [P(1) diff(P)];
   ulen = length(u)-1; % always even
   temp_x = [ones(ulen,1)*x(1); x(:); ones(ulen,1)*x(end)];
   x_smooth = conv(temp_x, u, 'valid');
   x_smooth = x_smooth(1+ulen/2:length(x)+ulen/2);
end

end

function p = normcdfBDP(x,mu,sigma)
% computes CDF of normal distribution using erfc()
z = (x-mu) ./ sigma;
p = NaN(size(z),class(z));

% Set edge case sigma=0
p(sigma==0 & x<mu) = 0;
p(sigma==0 & x>=mu) = 1;

% Normal cases
if isscalar(sigma)
   if sigma>0
      todo = true(size(z));
   else
      return;
   end
else
   todo = sigma>0;
end
z = z(todo);

% Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
% to produce accurate near-zero results for large negative x.
p(todo) = 0.5 * erfc(-z ./ sqrt(2));
end
