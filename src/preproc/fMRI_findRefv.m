function v_ref = fMRI_findRefv(v,k,verbose)
% v = nifti file
% k = searches only multiple of k timepoints
% verbose: 1 or 0
% author: soyoung choi USC


if isstruct(v)
    vls = v;
else isa(v,'char');
    if exist(v,'file')
        vls = load_untouch_nii_gz(v);
    end
end

n= floor(size(vls.img,4)/k);
smeasure = zeros(n,1);
tic
for i = 1: n
    t = i*k;
    s = motionEval(vls, t);
    smeasure(i) = mean(s);    
	if verbose
    		disp([num2str(t),': ',num2str(smeasure(i)),'; ',num2str(toc/60),' minutes'])
	end
end

[~,v_ref] = max(smeasure);
v_ref = v_ref*k;
end