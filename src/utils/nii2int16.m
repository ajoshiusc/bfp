
function nii2int16(in_nii, out_nii, negzero)

if ~exist('negzero','var')
    negzero = 0;
end

v=load_nii_BIG_Lab(in_nii);

v.hdr.dime.datatype=4; v.hdr.dime.bitpix = 16;

%amin=min(v.img(:)); v.img = v.img - amin; 
amax = max(v.img(:)); v.img = 1000*v.img/amax;

if negzero
    v.img(v.img<0)=0;
end

save_untouch_nii(v, out_nii);



