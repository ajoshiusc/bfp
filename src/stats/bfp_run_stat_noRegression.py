#%%
#config_file = '/Users/sychoi/Library/CloudStorage/Dropbox/FRIEND/Analysis/HBCD_NASCAR/V02_group3-2/config_BrainSync_atlas.ini'
import os
import sys
sys.stdout.flush()
config_file = os.environ['CONFIG_FILE']
#%%#%%
### Import the required librariesimport configparser
import sys
import os
import scipy.io as spio
import scipy as sp
import numpy as np
import configparser
import csv
from sklearn.linear_model import LinearRegression

### Import BrainSync libraries
config = configparser.ConfigParser()
config.read(config_file)
section = config.sections()
bfp_path = config.get('inputs','bfp_path')
sys.path.append(os.path.join(bfp_path, 'src/stats/') )
sys.path.append(os.path.join(str(bfp_path), 'src/BrainSync/')) 
from read_data_utils import load_bfp_dataT, read_demoCSV, write_text_timestamp,readConfig,load_bfp_dataT_dist2atlas
os.chdir(bfp_path)
cf = readConfig(config_file)
from brainsync import IDrefsub_BrainSync, groupBrainSync, generate_avgAtlas
from stats_utils import dist2atlas, sync2atlas, multiLinReg_corr
from stats_utils import randpairs_regression, multiLinReg_resid, LinReg_resid, multiLinReg_corr
from grayord_utils import vis_grayord_sigcorr, vis_grayord_sigpval
#%% 
if not os.path.isdir(cf.out_dir):
    os.makedirs(cf.out_dir)
log_fname = os.path.join(cf.out_dir, 'bfp_stat.log')
write_text_timestamp(log_fname, 'Config file used: ' + config_file +"\n All outputs will be written in: " + cf.out_dir )
# read demographic csv file
sub_ID, sub_fname, subAtlas_idx, reg_var, reg_cvar1, reg_cvar2 = read_demoCSV(cf.csv_fname,
                cf.data_dir,
                cf.file_ext,
                cf.colsubj,
                cf.colvar_exclude,
                cf.colvar_atlas,
                cf.colvar_main,
                cf.colvar_reg1,
                cf.colvar_reg2,
                cf.matcht,
                cf.lentime)
#%% makes file list for subjects
print('Identifying subjects for atlas creation and hypothesis testing...')
subTest_fname = []; subTest_IDs = []; subAtlas_fname = []; subAtlas_IDs = []
for ind in range(len(sub_ID)):
    sub = sub_ID[ind]
    fname = sub_fname[ind]
    if cf.stat_test == 'atlas-linear' or cf.stat_test == 'atlas-group':
        if int(subAtlas_idx[ind]) ==1:
            subAtlas_fname.append(fname)  
            subAtlas_IDs.append(sub)        
        else:
            subTest_fname.append(fname)
            subTest_IDs.append(sub)
    else:
        subTest_fname.append(fname)
        subTest_IDs.append(sub)


if cf.test_all == 'False':
    numT = len(subTest_IDs)
    write_text_timestamp(log_fname, "User Option: Only subjects not used for atlas creation will be used for hypothesis testing")            
else:
    numT = len(sub_ID)
    subTest_fname = sub_fname
    subTest_IDs = sub_ID
    write_text_timestamp(log_fname, "User Option: All subjects will be used for hypothesis testing")

count1=-1
subTest_varmain = np.zeros(numT); subTest_varc1 = np.zeros(numT); subTest_varc2 = np.zeros(numT)
for ind in range(len(sub_ID)):
    varmain = reg_var[ind]
    if len(reg_cvar1) != 0:
        varc1 = reg_cvar1[ind]
        varc2 = reg_cvar2[ind]
    if cf.test_all == 'False':
        if int(subAtlas_idx[ind]) !=1:  
            count1 +=1
    else:
        count1 +=1
    subTest_varmain[count1] = varmain
    
    if 'varc1' in locals():
        subTest_varc1[count1] = varc1
        subTest_varc2[count1] = varc2

del sub_ID, sub_fname, subAtlas_idx, reg_var, reg_cvar1, reg_cvar2, fname, sub, count1, ind, numT
if cf.stat_test == 'atlas-linear' or cf.stat_test == 'atlas-group':
    write_text_timestamp(log_fname, str(len(subAtlas_IDs)) + ' subjects will be used for atlas creation.' +
                         '\n '+ str(len(subTest_IDs)) + ' subjects will be used for hypothesis testing.')
#%%
# reads reference data and creates atlas by BrainSync algorithm
if cf.stat_test == 'atlas-linear' or cf.stat_test == 'atlas-group':
    if len(cf.atlas_fname) !=0:
        write_text_timestamp(log_fname, 'User Option: '+ cf.stat_test +
            '\n User defined atlas will be used ' + cf.atlas_fname)
        df = spio.loadmat(cf.atlas_fname)
        atlas_data = df['atlas_data']
        del df
    else:
        subAtlas_data, subAtlas_numT = load_bfp_dataT(subAtlas_fname, int(cf.lentime),cf.matcht)
        with open(cf.out_dir + "/subjects_atlas.csv", 'w') as csvfile:
            csv.writer(csvfile).writerows(zip(subAtlas_IDs, subAtlas_numT))
            del subAtlas_numT
        if cf.atlas_groupsync == 'True':
            write_text_timestamp(log_fname, 'User Option: '+ cf.stat_test +
                                 '\n Group BrainSync algorithm will be used for atlas creation')
            atlas_data, _, _, _ = groupBrainSync(subAtlas_data)
        else:
            write_text_timestamp(log_fname, 'User Option: '+ cf.stat_test +
                                 'Representative subject will be used for atlas creation')
            subRef_data, subRef_num = IDrefsub_BrainSync(subAtlas_data)
            write_text_timestamp(log_fname, 'Subject number ' + str(subAtlas_IDs[subRef_num]) + ' will be used for atlas creation')
            atlas_data = generate_avgAtlas(subRef_data, subAtlas_data)
        del subAtlas_data
        cf.atlas_fname = os.path.join(cf.out_dir + '/atlas.mat')
        spio.savemat(cf.atlas_fname, {'atlas_data': atlas_data})
        write_text_timestamp(log_fname, 'Atlas saved out as '+os.path.join(cf.out_dir + '/atlas.mat')) 

# #%% atlas-based linear regression
# if cf.stat_test == 'atlas-linear':
#     #load, sync and calculates geodesic distances
#     subTest_diff, numT = load_bfp_dataT_dist2atlas(subTest_fname, cf.atlas_fname, int(cf.lentime), cf.matcht)
#     spio.savemat(os.path.join(cf.out_dir + '/dist2atlas.mat'), {'subTest_diff': subTest_diff})
    
#     #write out csv file of subjects tested and variables used for testing
#     with open(cf.out_dir + "/subjects_testing.csv", 'w') as csvfile:
#         if cf.stat_test == 'atlas-linear' or cf.stat_test == 'atlas-group':
#             csv.writer(csvfile).writerows(zip(subTest_IDs, numT,subTest_varmain, subTest_varc1, subTest_varc2))
#         else:
#             csv.writer(csvfile).writerows(zip(subTest_IDs, subTest_varmain, subTest_varc1, subTest_varc2))
    
#     # computes correlation after controlling for two covariates
#     rval, pval, pval_fdr, msg = multiLinReg_corr(subTest_diff, subTest_varmain, subTest_varc1, subTest_varc2, float(cf.sig_alpha), 'linear')
#     spio.savemat(os.path.join(cf.out_dir + '/' + cf.outname + '_rval.mat'), {'rval': rval})
#     spio.savemat(os.path.join(cf.out_dir + '/' + cf.outname + '_pval.mat'), {'pval': pval})
#     spio.savemat(os.path.join(cf.out_dir + '/' + cf.outname + '_pval_fdr.mat'), {'pval_fdr': pval_fdr})
#     write_text_timestamp(log_fname, 'Done runnning atlas-based linear regression. ' + msg)
    
#     #%% visualization of results
#     vis_grayord_sigcorr(pval, rval, float(cf.sig_alpha), cf.outname, cf.out_dir, int(cf.smooth_iter), cf.save_figures, cf.bfp_path, cf.fsl_path)
#     vis_grayord_sigcorr(pval_fdr, rval, float(cf.sig_alpha), cf.outname + '_fdr', cf.out_dir, int(cf.smooth_iter), bool(cf.save_figures), cf.bfp_path, cf.fsl_path)
#     write_text_timestamp(log_fname, 'BFP regression analysis complete!'+
#           '\n How to interpret results:'+
#           '\n pvalue labeled surfaces have colorbar limits 0 to ' + str(cf.sig_alpha) + '; colorbar class is jet reverse.'+ 
#           '\n rvalue labeled surfaces have colorbar limits -0.5 to +0.5; colorbar class is jet. ' +
#           '\n positive rvalues indicate that higher ' + cf.colvar_main +' is associated with lower similarity to the atlas (higher geodesic distance) '+
#           'after controlling for '+ cf.colvar_reg1+' and '+cf.colvar_reg2 +'.'+
#           '\n negative rvalues indicate that higher '+ cf.colvar_main +' is associated with higher similarity to the atlas (lower geodesic distance) '+
#           'after controlling for '+ cf.colvar_reg1+' and '+cf.colvar_reg2 +'.' +
#           '\n adjusted pvalues corrected for multiple correction using FDR (alpha='+str(cf.sig_alpha)+')'+
#           '\n see https://doi.org/10.1016/j.neuroimage.2018.01.058 for further details')
    
# #%%



#%% PAIRWISE LINEAR REGRESSION TEST
# Do Linear regression on the covariates
# if cf.stat_test == 'pairwise-linear':
#     if cf.pw_fdr == 'True':
#         msg = 'adjusted pvalues corrected for multiple comparisons using FDR (alpha='+str(cf.sig_alpha)+')'
#     else:
#         msg = 'adjusted pvalues corrected for multiple comparisons using maxT permutation testing (permutations='+str(cf.pw_perm)+'; alpha='+str(cf.sig_alpha)+')'
#     #print(msg)
    
#     write_text_timestamp(
#         log_fname,
#         'Performing pair-wise linear regression.' +
#         '\n Number of pairs measured: ' + cf.pw_pairs +
#         '\n '+ msg)
# #    subTest_varc12 = sp.zeros((subTest_varc1.shape[0], 2))
# #    subTest_varc12[:, 0] = subTest_varc1
# #    subTest_varc12[:, 1] = subTest_varc2
# #    regr = LinearRegression()
# #    regr.fit(subTest_varc12, subTest_varmain)
# #    pre = regr.predict(subTest_varc12)
# #    subTest_varmain2 = subTest_varmain - pre    

# # Compute pairwise distance and perform regression
#     corr_pval_max, corr_pval_fdr = randpairs_regression(
#         bfp_path=cf.bfp_path,
#         sub_files=subTest_fname,
#         reg_var=subTest_varmain,
#         num_pairs=int(cf.pw_pairs),  # 19900,
#         nperm=int(cf.pw_perm),
#         len_time=int(cf.lentime),
#         num_proc=1,
#         pearson_fdr_test=bool(cf.pw_fdr))
    # # saves out results
    # spio.savemat(os.path.join(cf.out_dir + '/' + cf.outname + '_corr_pval_max.mat'), {'corr_pval_max': corr_pval_max})
    # spio.savemat(os.path.join(cf.out_dir + '/' + cf.outname + '_corr_pval_fdr.mat'), {'corr_pval_fdr': corr_pval_fdr})
    # #Visualization of the results
    # vis_grayord_sigpval(corr_pval_max, float(cf.sig_alpha), 
    #                     surf_name=cf.outname + '_max',
    #                     out_dir=cf.out_dir,
    #                     smooth_iter=int(cf.smooth_iter),
    #                     bfp_path=cf.bfp_path,
    #                     fsl_path=cf.fsl_path,
    #                     save_png=bool(cf.save_figures))
    
    # vis_grayord_sigpval(corr_pval_fdr, float(cf.sig_alpha), 
    #                     surf_name=cf.outname + 'fdr',
    #                     out_dir=cf.out_dir,
    #                     smooth_iter=int(cf.smooth_iter),
    #                     bfp_path=cf.bfp_path,
    #                     fsl_path=cf.fsl_path,
    #                     save_png=bool(cf.save_figures))
    
    # write_text_timestamp(log_fname, 'BFP regression analysis complete! '+
    #                      '\n pvalue labeled surfaces have colorbar limits 0 to ' + str(cf.sig_alpha) + '; colorbar class is jet reverse.'+ 
    #                      '\n significance indicates that brain connectivty is associated with '+ cf.colvar_main +
    #                      '\n see https://doi.org/10.1016/j.neuroimage.2018.01.058 for further details')
