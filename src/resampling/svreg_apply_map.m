% SVReg: Surface-Constrained Volumetric Registration
% Copyright (C) 2019 The Regents of the University of California and the University of Southern California
% Created by Anand A. Joshi, Chitresh Bhushan, David W. Shattuck, Richard M. Leahy
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; version 2.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
% USA.

function svreg_apply_map(map_file,data_file,out_file,target_file,smoothness,datatype,bitpix,interp_type)

if exist('datatype','var')
    if ~isempty(datatype)
        if ischar(datatype)
            datatype=str2double(datatype);
        end
    end
end

if exist('interp_type','var')
    if isempty(interp_type)
        interp_type='nearest';
    end
else
    interp_type='nearest';
end

if exist('bitpix','var')
    if ~isempty(bitpix)
        if ischar(bitpix)
            bitpix=str2double(bitpix);
        end
    end
end

if exist('smoothness','var')
    if ~isempty(smoothness)
        if ischar(smoothness)
            smoothness=str2double(smoothness);
        end
    else
        smoothness=1;
    end
end
vMap=load_nii_BIG_Lab(map_file);

if exist('smoothness','var')
    if smoothness == 0
        smoothness = [];
    end
    if ~isempty(smoothness)
        vMap.img(:,:,:,1)=smooth3(vMap.img(:,:,:,1),'gaussian',3*smoothness,smoothness);
        vMap.img(:,:,:,2)=smooth3(vMap.img(:,:,:,2),'gaussian',3*smoothness,smoothness);
        vMap.img(:,:,:,3)=smooth3(vMap.img(:,:,:,3),'gaussian',3*smoothness,smoothness);
    end
end

if ~contains(data_file,'.eig.')
    vsubjFA=load_nii_BIG_Lab(data_file);
    vAtlasFA=load_nii_BIG_Lab(target_file);vt=vAtlasFA;
    
    if length(size(vsubjFA.img))==3
        vAtlasFA.img=interp3(single(vsubjFA.img),vMap.img(:,:,:,2),vMap.img(:,:,:,1),vMap.img(:,:,:,3),interp_type);
    end
    
    % apply it to 4D data
    if length(size(vsubjFA.img))==4
        vAtlasFA.hdr.dime.dim(1) = 4;
        vAtlasFA.hdr.dime.dim(5) = size(vsubjFA.img,4);
        vAtlasFA.img = zeros(size(vAtlasFA.img,1),size(vAtlasFA.img,2),size(vAtlasFA.img,3),size(vsubjFA.img,4),'like',single([]));
        for j = 1:size(vsubjFA.img,4)
            vAtlasFA.img(:,:,:,j)=interp3(single(vsubjFA.img(:,:,:,j)),vMap.img(:,:,:,2),vMap.img(:,:,:,1),vMap.img(:,:,:,3),interp_type);
        end
        
    end
    
    vAtlasFA.img(isnan(vAtlasFA.img)) = 0;
    if exist('datatype','var')
        if ~isempty(datatype)
            vAtlasFA.hdr.dime.datatype=datatype;
            vAtlasFA.hdr.dime.bitpix=bitpix;
        else
            vAtlasFA.hdr.dime.datatype=vsubjFA.hdr.dime.datatype;
            vAtlasFA.hdr.dime.bitpix=vsubjFA.hdr.dime.bitpix;
        end
    else
        vAtlasFA.hdr.dime.datatype=vsubjFA.hdr.dime.datatype;
        vAtlasFA.hdr.dime.bitpix=vsubjFA.hdr.dime.bitpix;
    end
    if exist('smoothness','var')
        if ~isempty(smoothness)
            vAtlasFA.img=double(vAtlasFA.img).*double(vt.img>0);
        end
    end
    
    save_untouch_nii_gz(vAtlasFA,out_file);
else
    in_base=tempname;%
    %   subbasename=data_file(1:strfind(data_file,'.dwi.')-1);
    out_base=[in_base,'.atlas'];
    %  out_base_orig=[subbasename,'.atlas'];
    eig2nifti(data_file,in_base);
    warp_DTI(map_file, [in_base '.L1.nii.gz'], [in_base '.L2.nii.gz'], [in_base '.L3.nii.gz'], [in_base '.V1.nii.gz'], [in_base '.V2.nii.gz'], [in_base '.V3.nii.gz'], out_base);
    nifti2eig([out_base '.L1.nii.gz'], [out_base '.L2.nii.gz'], [out_base '.L3.nii.gz'], [out_base '.V1.nii.gz'], [out_base '.V2.nii.gz'], [out_base '.V3.nii.gz'], out_file);
    %     v=load_untouch_eig_gz(data_file);
    %     [interpL1, interpL2, interpL3, interpV1, interpV2, interpV3] = interp_dti_data(vMap.img, v.img(:,:,:,1), v.img(:,:,:,2), v.img(:,:,:,3), v.img(:,:,:,4:6), v.img(:,:,:,7:9), v.img(:,:,:,10:12), false);
    %     regMap=load_nii_BIG_Lab(map_file);
    %     [J1,J2,J3] = jacobian_nii(regMap, 'no_smooth');
    %     [W1,W2,W3] = PPD(interpV1, interpV2, interpV3, J1, J2, J3);
    %     generate_eig_file(interpL1, interpL2, interpL3, W1, W2, W3, out_file);
end
