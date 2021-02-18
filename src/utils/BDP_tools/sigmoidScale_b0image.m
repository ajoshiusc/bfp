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


function [img_out, cs] = sigmoidScale_b0image(img_in)
% Rescales b=0 image using a sigmoid kind of scaling, based on lookup table with spline fit. 
% Usage:
%   sigmoid_scale_b0image 
%   [img_out, cs] = sigmoid_scale_b0image(img_in)
%
% When no input argument is used, it shows a plot for comparison with other sigmoid scaling
% functions. 
%
% img_in - input b=0 image usually scaled between [0 C]; C>1
% img_out- output image. It would be scaled in range [0 1]. 
% cs     - The piecewise polynomial form of the cubic spline interpolant for later use with ppval.
%          See X, Y below for more details.
%

X = [0 0.1 0.2 0.3 0.4 0.5   0.6   0.7   0.8   0.9  1.0  1.1   1.2 ];
Y = [0 0.1 0.2 0.3 0.4 0.495 0.585 0.665 0.735 0.8  0.85 0.885 0.9105];

d = diff(tanh(1.2:0.1:6));
Y_rest = cumsum([Y(end) d*((1-Y(end))/sum(d))]);

X = [X 1.3:0.1:6     6.1 6.2 6.3 6.4 6.5];
Y = [Y Y_rest(2:end) 1   1   1   1   1  ];

cs = spline(X, [1 Y 0]);

if nargin==0  % to see nature of scaling
   x = 0:0.005:7;
   t = ppval(cs, x);
   
   figure;
   plot(x, t, 'r')
   hold on
   plot(x, erf(x), 'b')
   hold on
   plot(x, tanh(x), 'g')
   legend('sigmoid-scale-b0image(x)','erf(x)','tanh(x)', 'Location', 'SouthEast')
   grid on
   
else
   
   img_out = img_in;
   msk = img_in>0.4 & img_in<=6;
   img_out(msk) = ppval(cs, img_in(msk));
   img_out(img_in>6) = 1; % to avoid even a tiny problem at end
end
end

