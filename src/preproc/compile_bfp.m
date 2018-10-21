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

% Compile auxillary binaries
mcc -m -v nii2int16.m
mcc -m -v resample2surf.m
mcc -m -v generateSurfGOrdfMRI.m
mcc -m -v generateVolGOrdfMRI.m
mcc -m -v tNLMPDFGOrdfMRI.m
mcc -m -v generateGOrdSCT.m
mcc -m -v combineSurfVolGOrdfMRI.m

% Main BFP binary
copyfile('../scripts/bfp_linux.sh', [bfp_bin_dir,'/bfp.sh']);
movefile('bfp', bfp_bin_dir);

% Auxillary binaries
copyfile('../scripts/nii2int16_linux.sh', [bfp_bin_dir,'/nii2int16.sh']);
movefile('nii2int16', bfp_bin_dir);
copyfile('../scripts/resample2surf_linux.sh', [bfp_bin_dir,'/resample2surf.sh']);
movefile('resample2surf', bfp_bin_dir);
copyfile('../scripts/generateSurfGOrdfMRI_linux.sh', [bfp_bin_dir,'/generateSurfGOrdfMRI.sh']);
movefile('generateSurfGOrdfMRI', bfp_bin_dir);
copyfile('../scripts/generateVolGOrdfMRI_linux.sh', [bfp_bin_dir,'/generateVolGOrdfMRI.sh']);
movefile('generateVolGOrdfMRI', bfp_bin_dir);
copyfile('../scripts/tNLMPDFGOrdfMRI_linux.sh', [bfp_bin_dir,'/tNLMPDFGOrdfMRI.sh']);
movefile('tNLMPDFGOrdfMRI', bfp_bin_dir);
copyfile('../scripts/generateGOrdSCT_linux.sh', [bfp_bin_dir,'/generateGOrdSCT.sh']);
movefile('generateGOrdSCT', bfp_bin_dir);
copyfile('../scripts/combineSurfVolGOrdfMRI_linux.sh', [bfp_bin_dir,'/combineSurfVolGOrdfMRI.sh']);
movefile('combineSurfVolGOrdfMRI', bfp_bin_dir);

copyfile('../../supp_data', [bfp_bin_dir,'/supp_data']);
copyfile('bfp_version.txt',bfp_bin_dir);

% Pack the build
tar(sprintf('bfp_%s.tar.gz', bfp_version),bfp_bin_dir);
rmdir(bfp_bin_dir, 's');


