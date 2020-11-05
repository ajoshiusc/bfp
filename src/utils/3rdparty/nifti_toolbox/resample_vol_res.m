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

function vs=resample_vol_res(v,res_out,size_out,mthd,extrapval)


% This is a helper function for svreg_resample

if ~exist('mthd','var')
    mthd='linear';
end

size_in=size(v.img);
len_t=0;
if length(size_in)==4
    len_t=size_in(4);
    size_in=size_in(1:3);
end

res_in=v.hdr.dime.pixdim(2:4);

if ~isempty(size_out)
    res_out(1)=res_in(1)*size_in(1)/size_out(1);
    res_out(2)=res_in(2)*size_in(2)/size_out(2);
    res_out(3)=res_in(3)*size_in(3)/size_out(3);
end

size_mr(1)=size_in(1)*res_in(1);
size_mr(2)=size_in(2)*res_in(2);
size_mr(3)=size_in(3)*res_in(3);

xxin=res_in(1)/2:res_in(1):(size_mr(1)-res_in(1)/2);
xxout=res_out(1)/2:res_out(1):(size_mr(1)-res_out(1)/2);
yyin=res_in(2)/2:res_in(2):(size_mr(2)-res_in(2)/2);
yyout=res_out(2)/2:res_out(2):(size_mr(2)-res_out(2)/2);
zzin=res_in(3)/2:res_in(3):(size_mr(3)-res_in(3)/2);
zzout=res_out(3)/2:res_out(3):(size_mr(3)-res_out(3)/2);

[xxI,yyI,zzI]=ndgrid(double(xxin),double(yyin),double(zzin));
[xxO,yyO,zzO]=ndgrid(double(xxout),double(yyout),double(zzout));
size_out=size(xxO);
vs=v;vs.img=[];
if len_t>0
    for t=1:len_t
        vs.img(:,:,:,t)=interp3(yyI,xxI,zzI,double(v.img(:,:,:,t)),yyO,xxO,zzO,mthd,extrapval);
    end
else
    vs.img=interp3(yyI,xxI,zzI,double(v.img),yyO,xxO,zzO,mthd,extrapval);
end
vs.hdr.dime.pixdim(2:4)=res_out;
vs.hdr.dime.dim(2:4)=size_out;
