function [data_out, Tr_pp, Tr, bin_center] = histeqSmooth(data, arg2, arg3, arg4, arg5)
% Performs histogram equilization using parzen window estimates of histogram and pchip based
% mapping of CDFs. It should produce smoother output images and more continueous and smooth
% intensity transformation than matlab based histeq.
%
% Usage:
%    [data_out, Tr_pp, Tr, bin_center] = histeqSmooth(data, nbins)
%    [data_out, Tr_pp, Tr, bin_center] = histeqSmooth(data, hgram_target) % NOT RECOMMENDED
%    [data_out, Tr_pp, Tr, bin_center] = histeqSmooth(data, data_target, nbins)
%    [data_out, Tr_pp, Tr, bin_center] = histeqSmooth(data, data_target, nbins, win_width)
%    [data_out, Tr_pp, Tr, bin_center] = histeqSmooth(data, hgram_target, nbins, win_width, nthreads) 
%            
% In last usage, hgram_target is generated using histcParzen() with same nbins & win_width. It
% reduces computation time by not computing hgram_target repeatedly in an iterative algorithm. 
%
% Inputs:
%    data - Input data, matrix/vector; Must be scaled in range [0 1]
%    nbins - Number equidistant bins for equilization
%    hgram_target - A vector representing target distribution. It is implicitly assumed that
%                   bin_centers are linspace(0,1,length(hgram_target)).
%    data_target - Instead of specifying target histogram/distribution, a target data can be
%                  specified. In this case, data_out will have same intensity distribution as
%                  data_target. This is recommended method. 
%    win_width - Width of parzen window to be used for estimating distribution of all input data
%
% Outputs:
%    data_out - Histogram equalized output 
%    Tr_pp - Piecewise polynomial form of the intensity transformation, for later use with ppval.
%    [Tr, bin_center] - Output transformation map (Tr) is defined at each bin_center, such that 
%                       Tr_pp = pchip(bin_center, Tr);
%

nthreads = 6;
sz = size(data);
data = double(data(:));
if max(data)>1 || min(data)<0
   error('data must have values in range [0 1].')
end

if nargin==2
   if numel(arg2)==1
      nbins = arg2;
   elseif isvector(arg2)
      hgram_target = (arg2/sum(arg2))*numel(data);
      nbins = length(arg2);
   else
      error('arg2 must either an integer or a vector of histogram count.')
   end
   win_width = round(nbins*0.05); % 5% of nbins
   
elseif nargin==3 || nargin==4
   data_target = double(arg2(:));
   nbins = arg3;
   if nargin==4
      win_width = arg4;
   else
      win_width = round(nbins*0.12); % 12% of nbins
   end
   
   if max(data_target)>1 || min(data_target)<0
      error('data_target must have values in range [0 1].')
   end
   [~, ~, hgram_target, bin_center] = histcParzen(data_target, [0 1], nbins, win_width, nthreads);
   hgram_target = (hgram_target/sum(hgram_target))*numel(data);
   
elseif nargin==5
   hgram_target = (arg2/sum(arg2))*numel(data);
   nbins = arg3;
   win_width = arg4;
   nthreads = arg5;
else
   error('Incorrect number of inputs. see usage.')
end

[h1, b1, hgram_in, bin_center] = histcParzen(data, [0 1], nbins, win_width, nthreads);
hgram_in = (hgram_in/sum(hgram_in))*numel(data); % to match numel(data) perfectly

% Histogram equalization
if nargin==2 && numel(arg2)==1 % histeqSmooth(data, nbins)
   
   % generate ideal uniform histgram estimated from parzen window estimate
   parzen_kern = gen_parzen_kernal(win_width+1);   
   k = ones(1,nbins);
   k = (k/sum(k))*numel(data);
   hgram_target = conv(k, parzen_kern, 'full');
   
   cum_in = estimateInvertibleCDF(hgram_in(:), 0.5);
   cum_out = estimateInvertibleCDF(hgram_target(:), 0.5);
   Tr = pchip(cum_out(1:end-1), bin_center(1:end-1), cum_in); % last point is redundant in ideal case
   
elseif nargin==2 && isvector(arg2) % histeqSmooth(data, hgram_target)
   hgram_in = h1;
   bin_center = b1;
   hgram_in = (hgram_in/sum(hgram_in))*numel(data);
   
   cum_out = [0; estimateInvertibleCDF(hgram_in(:), 0.5)];
   cum_in = [0; estimateInvertibleCDF(hgram_target(:), 0.5)];
   Tr = pchip(cum_out, bin_center, cum_in);
   
else % histeqSmooth(data, data_target, nbins)
   cum_in = estimateInvertibleCDF(hgram_in(:), 0.5);
   cum_out = estimateInvertibleCDF(hgram_target(:), 0.5);
   [C,ia] = unique(cum_out);
   B = bin_center(ia);
   Tr = pchip(C, B, cum_in);
end

Tr(Tr<0) = 0;
Tr(Tr>1) = 1;
Tr_pp = pchip(bin_center, Tr);
data_out = ppval(Tr_pp, data);
data_out = reshape(data_out, sz);

end

function [parzen_kern, parzen_ind] = gen_parzen_kernal(window_sz)

parzen_ind = linspace(0,3,window_sz);
parzen_kern = zeros(size(parzen_ind));
parzen_kern(parzen_ind<=1) = 0.5* (parzen_ind(parzen_ind<=1).^2);
parzen_kern(parzen_ind>1 & parzen_ind<=2) = 0.5*(-2*(parzen_ind(parzen_ind>1 & parzen_ind<=2).^2) ...
   + 6*parzen_ind(parzen_ind>1 & parzen_ind<=2)-3);
parzen_kern(parzen_ind>2) = 0.5*((3-parzen_ind(parzen_ind>2)).^2);
parzen_kern = parzen_kern/sum(parzen_kern);

end

function c = estimateInvertibleCDF(p, tol)
% Invertible CDF estimate by replacing flat regions (p<=tol) in CDF by linear segment
if ~exist('tol', 'var')
   tol = 1e-4;
end

cdf_usual = cumsum(p(:));
c = cdf_usual;

msk = (p<=tol);
CC = bwconncomp(msk);
if CC.NumObjects >= 1
   for k = 1:CC.NumObjects
      ind_st = max(1, CC.PixelIdxList{k}(1)-1);
      ind_end = min(length(p), CC.PixelIdxList{k}(end)+1);
      l = ind_end-ind_st+1;
      c(ind_st:ind_end) = linspace(cdf_usual(ind_st), cdf_usual(ind_end),l);
   end
end

end


