function [I1out, low, high] = normalize_intensity(I1, percentile_range, I1_mask)
% Normalize images in range [0, 1]. The upper and lower cutoff intensities are decided by percentile
% as given in PERCENTILE_RANGE (0-100) inside the Mask (if specified)

percentile_range = percentile_range/100;
I1 = double(I1);

if ~exist('I1_mask','var')
   I1_mask = ones(size(I1));
end

I1_m = I1(I1_mask>0);

I1_m = sort(I1_m(:), 'ascend');
low = I1_m(max(floor(length(I1_m)*percentile_range(1)),1));
high = I1_m(ceil(length(I1_m)*percentile_range(2)));

if low==high
   I1out = ones(size(I1));
else
   I1out = double((I1-low)/(high-low));
   I1out(I1out<0) = 0;
   I1out(I1out>1) = 1;
end

end
