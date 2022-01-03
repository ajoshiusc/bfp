
function surfdata2gord(subbasename, left_data, right_data, GOrdSurfIndFile, GOrdFile)
% Tells the compiler that parallel processing is required

GordSize=96854;

[pth, sub] = fileparts(subbasename);

sub_left = readdfs([subbasename,'.left.mid.cortex.svreg.dfs']);
sub_right = readdfs([subbasename,'.right.mid.cortex.svreg.dfs']);

atlas_left = readdfs(fullfile(pth,'atlas.left.mid.cortex.svreg.dfs'));
atlas_right = readdfs(fullfile(pth,'atlas.right.mid.cortex.svreg.dfs'));

atlas_left = map_data_flatmap(sub_left,left_data,atlas_left);
atlas_right = map_data_flatmap(sub_right,right_data,atlas_right);

load(GOrdSurfIndFile,'ind_left','ind_right');

left_GO_SCT = atlas_left(ind_left,:);
right_GO_SCT = atlas_right(ind_right,:);

NDim = size(atlas_right,2);
lsz = size(left_GO_SCT,1);
rsz = size(right_GO_SCT,1);
data = nan(GordSize,NDim);
data(1:lsz,:) = left_GO_SCT;
data(1+lsz:lsz+rsz,:) = right_GO_SCT;

save(GOrdFile,'data');
