function generateVolGOrdfMRI(GOrdVolIndFile,subbasename,fmri,GOrdVolFile)

load(GOrdVolIndFile);

invmap=load_nii_BIG_Lab([subbasename,'.svreg.inv.map.nii.gz']);
invmapx=invmap.img(:,:,:,1);
invmapy=invmap.img(:,:,:,2);
invmapz=invmap.img(:,:,:,3);

ind=~isnan(bci_vol_ind);

t1=load_nii_BIG_Lab([subbasename,'.bfc.nii.gz']);
t1res=t1.hdr.dime.pixdim(2:4);
f=load_nii_BIG_Lab(fmri);
fres=f.hdr.dime.pixdim(2:4);

gr_x=invmapx(bci_vol_ind(ind));gr_x=gr_x*fres(1)/t1res(1);
gr_y=invmapy(bci_vol_ind(ind));gr_y=gr_y*fres(2)/t1res(2);
gr_z=invmapz(bci_vol_ind(ind));gr_z=gr_z*fres(3)/t1res(3);

Vol_GO_fmri=zeros(length(bci_vol_ind),size(f.img,3));
for t=1:size(f.img,3)
    Vol_GO_fmri(ind,t)=interp3(f.img,gr_y(ind),gr_x(ind),gr_z(ind));
end

save(GOrdVolFile,'Vol_GO_fmri');

