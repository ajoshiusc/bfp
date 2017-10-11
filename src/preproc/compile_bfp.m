function compile_bfp(bfp_version)
restoredefaultpath;
addpath(genpath('../../src'));

fid = fopen('bfp_version.txt', 'w');
fprintf(fid, bfp_version);
fclose(fid);

if exist(sprintf('bfp_%s', bfp_version),'dir')
    rmdir(sprintf('bfp_%s', bfp_version));
end

mkdir(sprintf('bfp_%s', bfp_version))
mcc -m -v bfp.m
copyfile('../scripts/bfp_linux.sh', sprintf('./bfp_%s/bfp.sh', bfp_version));
movefile('bfp', sprintf('bfp_%s', bfp_version));
copyfile('../../supp_data', sprintf('bfp_%s/supp_data', bfp_version));
copyfile('bfp_version.txt',sprintf('bfp_%s', bfp_version));



