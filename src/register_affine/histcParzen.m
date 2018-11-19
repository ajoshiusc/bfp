function [hist_count, bin_center, hist_count_total, bin_center_total] = histcParzen(x, x_range, nbins, win_width, nthreads)
% Computes parzen window estimate of the histogram of input data. The histogram bins are uniformly
% spaced in the specified x_range. Cubic B-spline window with variable width is used as the parzen
% window. This function is also interface to mex function histogram_parzen_variable_size_multithread_double.c
%
% Usage: 
%    [hist_count, bin_center, hist_count_total, bin_center_total] = histcParzen(x, x_range, nbins)
%    [hist_count, bin_center, hist_count_total, bin_center_total] = histcParzen(x, x_range, nbins, win_width)
%    [hist_count, bin_center, hist_count_total, bin_center_total] = histcParzen(x, x_range, nbins, win_width, nthreads)
%
% Inputs: 
%   x - Input intensity data, matrix or vector 
%   x_range - A vector of length two, denoting lower and upper level of binning range.
%   nbins - Number of bins
%   win_width - Width of parzen window. Default set to 8
%   nthreads - Number of parallel threads to use for computation. Default to 4.
%
% Outputs:
% hist_count - A vector of length nbins which has count in each bin. Note that due to use parzen
%              window, count in last bins may be off. Use the output hist_count_total for more
%              accurate estimate. 
% bin_center - A vector of length nbins which represents the center of each bin for hist_count. 
% hist_count_total, bin_center_total - These outputs have accurate estimate of counts near the edges
%                                      of the bins. However, the length of the vectors is more than
%                                      nbins and the range of binning is also slightly larger than
%                                      specified x_range (due to use of parzen window).
%                                      length(hist_count_total) = nbins + 2*floor(win_width/2)
%

if ~exist('win_width', 'var')
   win_width = 8;
end

if ~exist('nthreads', 'var')
   nthreads = 4;
end

% normalize image to [0 1]
x = (double(x(:))-x_range(1))/(x_range(2)-x_range(1));

% add extra bins on both sides to commpensate for width of parzen window
extra_bins = floor(win_width/2);
total_bins = nbins + 2*extra_bins;
temp = linspace(0, 1, nbins);
x_low = 0 - extra_bins*temp(2);
x_high = 1 + extra_bins*temp(2);
hist_count_total = histogram_parzen_variable_size_multithread_double(double(x), double(x_low), double(x_high), ...
   double(total_bins), double(win_width), double(nthreads));

% Deal with extra bins - throw them away
bin_center_total = (linspace(x_low, x_high, total_bins)*(x_range(2)-x_range(1))) + x_range(1);
bin_center = bin_center_total(1+extra_bins : end-extra_bins);
hist_count = hist_count_total(1+extra_bins : end-extra_bins);

% % Sum count from extra bins
% hist_count(1) = sum(hist_count_total(1:1+extra_bins));
% hist_count(end) = sum(hist_count_total(end-extra_bins:end));

% make col vector
bin_center = bin_center(:);
bin_center_total = bin_center_total(:);

end
