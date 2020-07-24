function s = motionEval(v, ref)
%   calculates structural similarity measure between reference volume and
%   all volumes in 4D time series data.
%   required input: v = fname or nifti file structure
%                   ref = reference volume, volume# or 'mid'
%   author: soyoung choi USC
if isstruct(v)
    vls = v;
else isa(v,'char');
    if exist(v,'file')
        vls = load_untouch_nii_gz(v);
    end
end

[x,y,z,nvol] = size(vls.img);
s= zeros(nvol,1);

if strcmpi(ref,'mid')
    k = round(nvol/2);
    rvls = vls;
    rvls.img = rvls.img(:,:,:,k);
elseif strcmpi(ref,'dssim')
    rvls = vls;
elseif isnumeric(ref)
    k = ref;
    rvls = vls;
    rvls.img = vls.img(:,:,:,k);
elseif isstruct(ref)
    rvls = ref;
elseif isfile(ref)
    rvls = load_untouch_nii_gz(ref);
end

datac = class(rvls.img);
if ~isa(datac,class(vls.img))
    rvls.img = double(rvls.img);
    vls.img = double(vls.img);
end

for i = 1:nvol
    if strcmpi(ref,'dssim')
        data2 = vls.img(round(x/4):round(3*x/4),round(y/4):round(3*y/4),round(z/4):round(3*z/4),i);
        if i == nvol
            data1 = vls.img(round(x/4):round(3*x/4),round(y/4):round(3*y/4),round(z/4):round(3*z/4),i-1);
        else
            data1 = vls.img(round(x/4):round(3*x/4),round(y/4):round(3*y/4),round(z/4):round(3*z/4),i+1);
        end
    else
        data1 = rvls.img(round(x/4):round(3*x/4),round(y/4):round(3*y/4),round(z/4):round(3*z/4));
        data2 = vls.img(round(x/4):round(3*x/4),round(y/4):round(3*y/4),round(z/4):round(3*z/4),i);
    end
    s(i,1) = ssim(data1,data2);
end
end