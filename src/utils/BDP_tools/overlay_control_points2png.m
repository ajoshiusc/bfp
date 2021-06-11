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


function overlay_control_points2png(epiImgR_plot, O_trans, spacing_plot, fname)
% Plot controls points slice wise over image

x = O_trans(:,:,:,1)+1;
y = O_trans(:,:,:,2)+1;

x_max = max(x(:))+4;
x_min = min(x(:))-4;
y_max = max(y(:))+4;
y_min = min(y(:))-4;

h = figure('Units','pixels','Position',[10 10 1000 1000]);
for n = 1:size(O_trans, 3)
   imagesc(epiImgR_plot(:,:,(n-1)*spacing_plot(3)+1), [0 1]);
   colormap(gray)
   hold on
   
   x = O_trans(:,:,n,1)+1;
   y = O_trans(:,:,n,2)+1;

   ylim([x_min x_max]);
   xlim([y_min y_max]);
   set(gca,'Color',[0 0 0.7]);
   axis equal 
   
   plot(y(:),x(:),'.', 'MarkerEdgeColor','r', 'MarkerSize',8);
   set(gcf, 'InvertHardCopy', 'off');
   saveas(gcf, [fname '.' num2str(n) '.png'])
   clf
   hold off
end
close(h)

end
