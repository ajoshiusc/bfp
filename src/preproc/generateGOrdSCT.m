
function generateGOrdSCT(subbasename, GOrdSurfIndFile)
% Tells the compiler that parallel processing is required
%#function gcp

NLevels=5;
GordSize=96854;
GOrdFile=[subbasename,'.SCT.GOrd.mat'];

%% process left hemisphere
sub_surf = [subbasename,'.left.mid.cortex.svreg.dfs'];
subdir = fileparts(subbasename);
atlas_surf = fullfile(subdir,'/atlas.left.mid.cortex.svreg.dfs');
thickness_surf = [subbasename, '.pvc-thickness_0-6mm.left.mid.cortex.dfs'];
SCT_outfile_left = [subbasename,'.left.SCT.mat'];
if ~exist(SCT_outfile_left,'file')
    cortical_SCT_measures(sub_surf,atlas_surf,NLevels,thickness_surf,SCT_outfile_left);
end

%% process right hemisphere
sub_surf = [subbasename,'.right.mid.cortex.svreg.dfs'];
subdir = fileparts(subbasename);
atlas_surf = fullfile(subdir,'/atlas.right.mid.cortex.svreg.dfs');
thickness_surf = [subbasename, '.pvc-thickness_0-6mm.right.mid.cortex.dfs'];
SCT_outfile_right = [subbasename,'.right.SCT.mat'];
if ~exist(SCT_outfile_right,'file')
    cortical_SCT_measures(sub_surf, atlas_surf, NLevels, thickness_surf, SCT_outfile_right);
end

%%Load the GOrd and map to atlas

load(GOrdSurfIndFile,'ind_left','ind_right');
left_sct = load(SCT_outfile_left);
right_sct = load(SCT_outfile_right);

left_GO_SCT=left_sct.SCTatlas(ind_left,:);
right_GO_SCT=right_sct.SCTatlas(ind_right,:);

NDim = size(right_sct.SCTatlas,2);
lsz=size(left_GO_SCT,1);
rsz=size(right_GO_SCT,1);
SCT_GO=zeros(GordSize,NDim);
SCT_GO(1:lsz,:)=left_GO_SCT;
SCT_GO(1+lsz:lsz+rsz,:)=right_GO_SCT;

save(GOrdFile,'SCT_GO');
