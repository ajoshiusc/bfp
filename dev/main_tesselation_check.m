

s=readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kright.dfs');
[FaceNormal, FaceArea, FaceCenter, VertexNormal, VertexArea, SuspectFace, NumFacesEachVertex, Duplicated_Faces, Not_Twice_Faces] = tessellation_stats(s,1);


s=readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kleft.dfs');
[FaceNormal, FaceArea, FaceCenter, VertexNormal, VertexArea, SuspectFace, NumFacesEachVertex, Duplicated_Faces, Not_Twice_Faces] = tessellation_stats(s,1);
