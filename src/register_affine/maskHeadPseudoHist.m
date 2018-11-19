function [mask, thresh] = maskHeadPseudoHist(data_file, output_mask_file, if_aggresive)
% This function tries to detect head (including skull and skin) just based on
% intensity in a 3D volume. This function assumes that the head (or area of interest) and
% background have significant difference in their intensity values.  
%
% output_mask_file - (optional) file name of output mask
%
% Always check output - may NOT make any sense at all!
%

nbins = 300; 
parzen_width = 40;
nthreads = 6;
log_sz = 21; % should be < parzen_width
log_sgma = 1;
log_thresh = 400; 
valley_slp_thresh = 5;

% check input options
if nargin==1
   if_aggresive = false;
elseif nargin==2 && ~ischar(output_mask_file)
   if_aggresive = output_mask_file;
   clear output_mask_file
elseif nargin==2
   if_aggresive = false;
end

% load input image and some sanity test
if ischar(data_file)
   data = load_untouch_nii_gz(data_file, true);
else
   if ~isfield(data_file,'untouch') || data_file.untouch ~= 1
      error('Please use ''load_untouch_nii.m'' to load the file.');
   end
   data = data_file;
   clear data_file
end

if ndims(data.img)==4
   fprintf('4D volume found. Only first 3D-volume will be used for masking.');
   data.img = data.img(:,:,:,1);
   data.hdr.dime.dim(5) = 1;
   data.hdr.dime.dim(1) = 3;
elseif ndims(data.img)>4
   error('data_file must be 3D volume.')
end

% check if the input is a mask itself
unq_intensities = unique(data.img(:));
if length(unq_intensities)<=2
   mask = data;

   if length(unq_intensities)==2
      threshold = min(unq_intensities);
   elseif length(unq_intensities)==1
      threshold = unq_intensities-1;
   end
   mask.img = mask.img>threshold;
   
   % make BrainSuite compatible
   mask.img = uint8(mask.img>0)*255;
   
   if exist('output_mask_file', 'var')
      mask.hdr.dime.scl_slope = 0;
      save_untouch_nii_gz(mask, output_mask_file, 2);
   end
   
   return;
end

[dataN, l, h] = normalize_intensity(data.img, [0 98]); % lower percentile should be 0 for robustness

% get thrsholds based histogram
[~,~, hist_data, bin_cntr] = histcParzen(dataN(:), [0 1], nbins, parzen_width, nthreads);

% -ve 2nd derivative = Peak
% +ve 2nd derivative = Valley
gh = gradient(hist_data);
[gh0cross, gh0cross_slp] = zeroCrossing(gh, bin_cntr);
gh0val = interp1(bin_cntr, hist_data, gh0cross);

[~, max_ind] = max(gh0val);
max_peak = gh0cross(max_ind); % histogram peak should be approx background intensity

if if_aggresive % use 1st valley after the max_peak
   ind = find((gh0cross_slp>valley_slp_thresh) & (gh0cross>=max_peak));
   if isempty(ind)
      threshold = max_peak;
   else
      threshold = gh0cross(ind(1)); % 1st large zero crossing with +ve slope
   end
   
else % use zero corssing of 2nd derivative of histogram
   lg_kr = log_kernel([1 log_sz], log_sgma);
   dd_hist = conv(hist_data, lg_kr,'same');
   [bin0cross, bin0cross_slp] = zeroCrossing(dd_hist, bin_cntr);
   ind = find(bin0cross_slp>log_thresh);
   if isempty(ind)
      threshold = max_peak;
   else
      threshold = bin0cross(ind(1)); % 1st large zero crossing with +ve slope
   end
end

% sanity check for detected threshold 
temp = dataN>threshold;
vol_ratio = 100*sum(temp(:))/numel(temp);
if vol_ratio<5 || threshold<0 % brain should be more than 5% of total volume 
   threshold = max_peak;
end

mask_data = dataN>threshold;
mask_data = imfill(mask_data, 'holes'); % close holes 


% cleanup mask - remove isolated components

% Structuring elements
se1 = strel_sphere(1);
se2 = strel_sphere(2);
se3 = strel_sphere(3);
se4 = strel_sphere(4);
data_size = size(mask_data);

% padd zeros to avoid weird behaviour at boundaries
pad_size = 4;
mask_data_pad = padarray(mask_data, [1 1 1]*pad_size)>0;

% remove isolated pixels
msk_tmp = imerode(mask_data_pad, se2);
msk_tmp = largest_connected_component_MASKHEADPSEUDOHIST(msk_tmp, 6);
msk_tmp = imdilate(msk_tmp, se3);
msk_tmp = imdilate(msk_tmp, se2);
msk_tmp = msk_tmp & mask_data_pad;
msk_tmp = largest_connected_component_MASKHEADPSEUDOHIST(msk_tmp, 6);
msk_tmp = imfill(msk_tmp>0, 'holes');

mask_data = msk_tmp(pad_size+1:data_size(1)+pad_size, pad_size+1:data_size(2)+pad_size, pad_size+1:data_size(3)+pad_size)>0;
clear msk_tmp mask_data_pad

% update threshold in terms of input data
thresh = l + threshold*(h-l);


% make BrainSuite compatible
mask = data;
mask.img = uint8(mask_data)*255;

if exist('output_mask_file', 'var')
   mask.hdr.dime.scl_slope = 0;   
   save_untouch_nii_gz(mask, output_mask_file, 2);
end

end


% Find the biggest connected component
function out = largest_connected_component_MASKHEADPSEUDOHIST(mask, conn)

cc = bwconncomp(mask, conn);
if cc.NumObjects>1
   cc_size = [];
   for k = 1:cc.NumObjects
      cc_size(k) = length(cc.PixelIdxList{k});
   end
   [~,IX] = sort(cc_size,'descend');
   out = false(size(mask));
   out(cc.PixelIdxList{IX(1)}) = true;
else
   out = mask;
end
end


function [bin0cross, bin0cross_slp, bin0kiss] = zeroCrossing(x, bin)
% find zero crossing of function x defined at bin

n = length(x);
ind0_kiss = [];

% First handle exact zeros
ind0 = find(x==0);
if ~isempty(ind0)
   tempind = ind0;
   if tempind(1)==1
      tempind(1) = [];
      ind0_kiss = [ind0_kiss 1];
   end
   if tempind(end)==n
      tempind(end) = [];
      ind0_kiss = [ind0_kiss n];
   end
   
   % following misses contiguous zero values of x, but thats OK here
   px = x(tempind+1) .* x(tempind-1);
   ind0_kiss = [ind0_kiss tempind(px>0)];
end

% look for zero crossings
ap = x(1:end-1) .* x(2:end); % adjacent points should have different signs
ind1 = find(ap<0);

frc = abs(x(ind1) ./ (x(ind1+1) - x(ind1)));
bin0cross = bin(ind1) + (frc .* (bin(ind1+1) - bin(ind1))); % actual location of zero crossing

% add exact zeros & sort
bin0cross = sort([bin0cross(:); vect(bin(ind0))], 'ascend');
bin0kiss = sort(vect(bin(ind0_kiss)), 'ascend');

% compute gradients 
dx = gradient(x, bin(2)-bin(1));
bin0cross_slp = interp1(bin, dx, bin0cross);

end

function t = prcentileIntensity(d, prct)
d = sort(d(:), 'ascend');
ind = ceil(prct/100*length(prct));
t = d(ind);
end

function img = stretch_intensity(img, msk)
if nargin==1
   msk = true(size(img));
end

l = min(img(msk));
h = max(img(msk));
img = (img-l)/(h-l);
end

function h = log_kernel(p2, p3)
% Laplacian of Gaussian
% h = log_kernel(size, sigma)

% first calculate Gaussian
siz   = (p2-1)/2;
std2   = p3^2;

[x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
arg   = -(x.*x + y.*y)/(2*std2);

h     = exp(arg);
h(h<eps*max(h(:))) = 0;

sumh = sum(h(:));
if sumh ~= 0,
   h  = h/sumh;
end;
% now calculate Laplacian
h1 = h.*(x.*x + y.*y - 2*std2)/(std2^2);
h     = h1 - sum(h1(:))/prod(p2); % make the filter sum to zero

end
