function cortical_SCT_measures(sub_surf,atlas_surf,NLevels,thickness_surf, SCT_outfile)

if ~exist('NLevels','var')
    NLevels=5;
end

if ~exist('thickness_surf','var')
    thickness_surf=sub_surf;
end

s=readdfs(sub_surf);
s=smooth_cortex_fast(s,.1,1000);
% compute principal curvatures, mean curvature, gaussian curvature
tic
[Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=patchcurvature(s,1);
toc
tmp1 = max (Lambda1, Lambda2);
tmp2 = min (Lambda1, Lambda2);
Lambda1 = tmp1;
Lambda2 = tmp2;

p=prctile(Lambda1,[1,99]);
Lambda1=max(p(1),min(Lambda1,p(2)));
p=prctile(Lambda2,[1,99]);
Lambda2=max(p(1),min(Lambda2,p(2)));
p=prctile(Cmean,[1,99]);
Cmean=max(p(1),min(Cmean,p(2)));
p=prctile(Cgaussian,[1,99]);
Cgaussian=max(p(1),min(Cgaussian,p(2)));

a=bipolar;
figure;
patch('faces',s.faces,'vertices',s.vertices,'facevertexcdata',-Cmean,'edgecolor','none','facecolor','interp');
caxis([-.25,.25]); axis equal;camlight;material dull;view(-90,0);camlight;colormap(a);colorbar;axis off;
figure;
patch('faces',s.faces,'vertices',s.vertices,'facevertexcdata',Cgaussian,'edgecolor','none','facecolor','interp');
caxis([-.0015,.0015]); axis equal;camlight;material dull;view(-90,0);camlight;colormap(a);colorbar;axis off;

% compute shape index and curvedness
S=-(2/pi)*atan((Lambda2+Lambda1)./(Lambda2-Lambda1));
C=0.5*(Lambda2.^2+Lambda1.^2).^0.5 ;

a=bipolar;
figure;
patch('faces',s.faces,'vertices',s.vertices,'facevertexcdata',S,'edgecolor','none','facecolor','interp');
caxis([-1,1]); axis equal;camlight;material dull;view(-90,0);camlight;colormap(a);colorbar;axis off;
figure;
patch('faces',s.faces,'vertices',s.vertices,'facevertexcdata',C,'edgecolor','none','facecolor','interp');
caxis([0,.5]); axis equal;camlight;material dull;view(-90,0);camlight;colorbar;axis off;

% shape index at multiple scales
Sm(:,1)=S;
for jj=2:NLevels
    Sm(:,jj)=smooth_surf_function(s,Sm(:,jj-1));
    figure;
    patch('faces',s.faces,'vertices',s.vertices,'facevertexcdata',Sm(:,jj),'edgecolor','none','facecolor','interp');
    caxis([-1,1]); axis equal;camlight;material dull;view(-90,0);camlight;colormap(a);colorbar;axis off;
end
% curvedness at multiple scales
Cm(:,1)=C;
for jj=2:NLevels
    Cm(:,jj)=smooth_surf_function(s,Cm(:,jj-1));
    figure;
    patch('faces',s.faces,'vertices',s.vertices,'facevertexcdata',Cm(:,jj),'edgecolor','none','facecolor','interp');
    caxis([0,.5]); axis equal;camlight;material dull;view(-90,0);camlight;colormap(a);colorbar;axis off;
end
% thickness at multiple scales

th_surf = readdfs(thickness_surf);
Tm(:,1) = th_surf.attributes;
for jj=2:NLevels
    Tm(:,jj)=smooth_surf_function(s,Tm(:,jj-1));
    figure;
    patch('faces',s.faces,'vertices',s.vertices,'facevertexcdata',Tm(:,jj),'edgecolor','none','facecolor','interp');
    caxis([0,4]); axis equal;camlight;material dull;view(-90,0);camlight;colormap(a);colorbar;axis off;
end

% Sm is a 2D matrix: (number of mesh nodes X number of scales)
SCT=[Sm,Cm,Tm];
% SCT is a 2D matrix: (number of mesh nodes X thrice the number of scales)

% warp SCT to atlas coordinate space
at=readdfs(atlas_surf);
SCTatlas=zeros(size(at.vertices,1),size(SCT,2));
for jj=1:size(SCT,2)
    SCTatlas(:,jj)=map_data_flatmap(s,SCT(:,jj),at);
    fprintf('mapping %d/%d SCT measure to atlas\n',jj,size(SCT,2));
end

save(SCT_outfile,'SCT','SCTatlas');

