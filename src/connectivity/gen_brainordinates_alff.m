function gen_brainordinates_alff(BSTDir, anatbasename, fmribasename,measure)

% Measure could be 'ALFF', 'ALFF_Z', or 'fALFF_Z'
%bfp_out,subid,'func',[subid,'_',sessionid
sub_mri = [anatbasename, '.bfc.nii.gz'];
subfmri = [fmribasename,'_',measure,'2standard.nii.gz'];
if ~isfile(subfmri)
    fprintf('File doesn''t exist:%s\n',subfmri);
    return;
end
map_file = [anatbasename,'.svreg.inv.map.nii.gz'];
mappedfmri = [fmribasename,'_bold_',measure,'2standard.mapped.nii.gz'];
bordfmri = [fmribasename,'_bold_',measure,'.BOrd.mat'];

if exist(bordfmri,'file')
    return;
end

%fullfile(funcDir,sprintf('%s_%s_bold.nii.gz',subid,sessionid{ind})),'file')

target = fullfile(BSTDir, 'svreg','BCI-DNI_brain_atlas','BCI-DNI_brain.pvc.frac.nii.gz');
%app_map_exe = fullfile(BSTDir,'svreg','bin','svreg_apply_map.sh');
%target_mask='/home/ajoshi/BrainSuite19b/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.mask.nii.gz';

resampled_target = [tempname(),'.nii.gz'];
resampled_map =  [tempname(),'.nii.gz'];
% resample mask to 3mm
%svreg_resample_exe = fullfile(BSTDir,'svreg','bin','svreg_resample.sh');
%cmd = sprintf('%s %s %s', svreg_resample_exe, target_mask, resampled_mask);
%unix(cmd);


svreg_resample(target, resampled_target,'-res','3','3','3')
svreg_resample(map_file, resampled_map,'-res','3','3','3')

v=load_nii_BIG_Lab(sub_mri);
mri_res=v.hdr.dime.pixdim(2:4);clear v;
map=load_nii_BIG_Lab(resampled_map);
map.img(:,:,:,1)=map.img(:,:,:,1)*mri_res(:,1)/3;
map.img(:,:,:,2)=map.img(:,:,:,2)*mri_res(:,2)/3;
map.img(:,:,:,3)=map.img(:,:,:,3)*mri_res(:,3)/3;
save_untouch_nii_gz(map,resampled_map);

svreg_apply_map(resampled_map, subfmri, mappedfmri, resampled_target);

vv=load_nii_BIG_Lab(resampled_target);
ind = find(vv.img(:)>0);

data=load_nii_BIG_Lab(mappedfmri);
numt = size(data.img,4);
dtseries = zeros(length(ind),numt);

for j=1:numt
    tmp=data.img(:,:,:,j);
    dtseries(:,j)=tmp(ind);
end

save(bordfmri,'dtseries');

delete(resampled_target);
delete(resampled_map);
