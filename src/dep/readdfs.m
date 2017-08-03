% SVReg: Surface-Constrained Volumetric Registration
% Copyright (C) 2017 The Regents of the University of California and the University of Southern California
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



function [NFV,hdr]=readdfs(s)
% updated 21-April-2010
% WRITEDFS2 Writes a Duff Surface file (dfs).
%   WRITEDFS2(FILENAME,NFV) writes the file specified by FILENAME string.
%
% DFS has the following structure:
% expected fields:
%   NFV.faces        : the face data,
%   NFV.vertices     : the vertex data
% optional fields:
%   NFV.normals      : surface normals
%   NFV.uv           : the uv coordinates,
%   NFV.vcolor       : the vertex color data (r,g,b) in ([0-1],[0-1],[0-1])
%   NFV.labels       : vertex label (int16)
%   NFV.attributes   : vertex attribute (float)
%
% Author : David Shattuck (shattuck@loni.ucla.edu)
% updated 21-April-2010

% Header format
% [000-011]	char headerType[12]; // should be DFS_BE v2.0\0 on big-endian machines, DFS_LEv1.0\0 on little-endian
% [012-015] int32 hdrsize;			// Size of complete header (i.e., offset of first data element)
% [016-010] int32 mdoffset;			// Start of metadata.
% [020-023] int32 pdoffset;			// Start of patient data header.
% [024-027] int32 nTriangles;		// Number of triangles
% [028-031] int32 nVertices;		// Number of vertices
% [032-035] int32 nStrips;			// Number of triangle strips (deprecated)
% [036-039] int32 stripSize;		// size of strip data  (deprecated)
% [040-043] int32 normals;			// 4	Int32	<normals>	Start of vertex normal data (0 if not in file)
% [044-047] int32 uvoffset;			// Start of surface parameterization data (0 if not in file)
% [048-051] int32 vcoffset;			// vertex color
% [052-055] int32 labelOffset;	// vertex labels
% [056-059] int32 vertexAttributes; // vertex attributes (float32 array of length NV)
% [060-183] uint8 pad2[4 + 15*8]; // formerly 4x4 matrix, affine transformation to world coordinates, now used to add new fields
if strcmp(s(end-2:end),'.gz')
    gunzip(s);
    s=s(1:end-3);s1=s;
end
fid=fopen(s,'rb','ieee-le');
if (fid<0)
    error('Unable to read:%s',s);
end
hdr.magic=char(fread(fid,12,'char'));
%disp(hdr.magic');
hdr.hdrsize=fread(fid,1,'int32');
hdr.mdoffset=fread(fid,1,'int32');
hdr.pdoffset=fread(fid,1,'int32');
hdr.nTriangles=fread(fid,1,'int32');
hdr.nVertices=fread(fid,1,'int32');
hdr.nStrips=fread(fid,1,'int32');
hdr.stripSize=fread(fid,1,'int32');
hdr.normals=fread(fid,1,'int32');
hdr.uvStart=fread(fid,1,'int32');
hdr.vcoffset=fread(fid,1,'int32');
hdr.labelOffset=fread(fid,1,'int32');
hdr.vertexAttributes=fread(fid,1,'int32');

fseek(fid,hdr.hdrsize,-1);
NFV.faces = fread(fid,[3 hdr.nTriangles],'int32')+1;
NFV.vertices=fread(fid,[3 hdr.nVertices],'float32');
NFV.faces=NFV.faces';
NFV.vertices=NFV.vertices';
if (hdr.normals>0)
    %display('reading vertex normals.');
    fseek(fid,hdr.normals,-1);
    NFV.normals = fread(fid,[3 hdr.nVertices],'float32')';
end;
if (hdr.vcoffset>0)
    %display('reading vertex colors.');
    fseek(fid,hdr.vcoffset,-1);
    NFV.vcolor = fread(fid,[3 hdr.nVertices],'float32')';
end;
if (hdr.uvStart>0)
    %display('reading uv coordinates.');
    fseek(fid,hdr.uvStart,-1);
    uv = fread(fid,[2 hdr.nVertices],'float32');
    NFV.u = uv(1,:);
    NFV.v = uv(2,:);
end;
if (hdr.labelOffset>0)
    %display('reading vertex labels.');
    fseek(fid,hdr.labelOffset,-1);
    NFV.labels = fread(fid,[hdr.nVertices],'uint16');
end;
if (hdr.vertexAttributes>0)
    %display('reading vertex attributes.');
    fseek(fid,hdr.vertexAttributes,-1);
    NFV.attributes = fread(fid,[hdr.nVertices],'float32');
end;
NFV.name=s;
fclose(fid);
if exist('s1','var')
    delete(s1);
end
