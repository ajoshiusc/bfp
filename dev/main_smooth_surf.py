import os
import sys
sys.path.append('/ImagePTE1/ajoshi/code_farm/bfp/src/stats')
from surfproc import view_patch_vtk, smooth_patch
from dfsio import readdfs,writedfs
bfp_path='/ImagePTE1/ajoshi/code_farm/bfp'
lsurf=''
smooth_iter=1000
lsurf = readdfs(os.path.join(bfp_path, 'supp_data/bci32kleft.dfs'))
lsurf = smooth_patch(lsurf, iterations=int(smooth_iter))
writedfs(os.path.join(bfp_path, 'supp_data/bci32kleft_smooth.dfs'),lsurf)

lsurf = readdfs(os.path.join(bfp_path, 'supp_data/bci32kright.dfs'))
lsurf = smooth_patch(lsurf, iterations=int(smooth_iter))
writedfs(os.path.join(bfp_path, 'supp_data/bci32kright_smooth.dfs'),lsurf)

