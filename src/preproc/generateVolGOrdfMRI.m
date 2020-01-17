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

gr_x=invmapx(bci_vol_ind(ind));gr_x=gr_x*t1res(1)/fres(1);
gr_y=invmapy(bci_vol_ind(ind));gr_y=gr_y*t1res(2)/fres(2);
gr_z=invmapz(bci_vol_ind(ind));gr_z=gr_z*t1res(3)/fres(3);

Vol_GO_fmri=zeros(length(bci_vol_ind),size(f.img,4));
for t=1:size(f.img,4)
    Vol_GO_fmri(ind,t)=interp3(f.img(:,:,:,t),gr_y,gr_x,gr_z);
end

save(GOrdVolFile,'Vol_GO_fmri');

