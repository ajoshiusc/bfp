function compile_bfp(bfp_version)
restoredefaultpath;
addpath(genpath('../../src'));

bfp_bin_dir=sprintf('bfp_%s', bfp_version);
fid = fopen('bfp_version.txt', 'w');
fprintf(fid, bfp_version);
fclose(fid);

if exist(bfp_bin_dir,'dir')
    rmdir(bfp_bin_dir,'s');
end

mkdir(bfp_bin_dir);
mcc -m -v bfp.m
mcc -m -v nii2int16.m
mcc -m -v resample2surf.m
mcc -m -v generateGOrdSurfFile.m

copyfile('../scripts/bfp_linux.sh', [bfp_bin_dir,'/bfp.sh']);
movefile('bfp', bfp_bin_dir);
copyfile('../scripts/nii2int16_linux.sh', [bfp_bin_dir,'/nii2int16.sh']);
movefile('nii2int16', bfp_bin_dir);

copyfile('../../supp_data', [bfp_bin_dir,'/supp_data']);
copyfile('bfp_version.txt',bfp_bin_dir);

tar(sprintf('bfp_%s.tar.gz', bfp_version),bfp_bin_dir);
rmdir(bfp_bin_dir, 's');


