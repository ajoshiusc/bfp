function C = convnPadded(A, B)
% Applies replicate padding before convolution. Output is same size as A.

padsize = (size(B)-1)/2;
padsize_pre = floor(padsize);
padsize_post = ceil(padsize);

A_padded = padarray(A, padsize_pre, 'replicate', 'pre');
A_padded = padarray(A_padded, padsize_post, 'replicate', 'post');

C = convn(A_padded, B, 'valid');

end
