def bfp_lincorr(config_file):
    #%%
    config_file = '/home/sychoi/Documents/MATLAB/bfp/src/stats/sample_config_stats.ini'
    #%%#%%
    ### Import the required librariesimport configparser
    import sys
    import os
    import scipy.io as spio
    import scipy as sp
    import numpy as np
    import configparser

    ### get BFP directory from config file
    config = configparser.ConfigParser()
    config.read(config_file)
    section = config.sections()
    bfp_path = config.get('inputs','bfp_path')
    sys.path.append(os.path.join(bfp_path, 'src/stats/') )
    from dev_config import readConfig
    cf = readConfig(config_file)

    ### Import BrainSync libraries
    sys.path.append(os.path.join(str(cf.bfp_path), 'src/BrainSync/')) 
    from brainsync import IDrefsub_BrainSync, groupBrainSync, generate_avgAtlas
    sys.path.append(cf.bfp_path + "/src/stats/") 
    from stats_utils import dist2atlas, load_bfp_data, read_demoCSV, sync2atlas, multiLinReg_corr
    from grayord_utils import vis_grayord_sigcorr
    #%% 
    os.chdir(bfp_path)
    print(cf.out_dir + ": writing output directory")
    if not os.path.isdir(cf.out_dir):
        os.makedirs(cf.out_dir)
    # read demographic csv file
    sub_ID, sub_fname, subAtlas_idx, reg_var, reg_cvar1, reg_cvar2 = read_demoCSV(cf.csv_fname,
                    cf.data_dir,
                    cf.file_ext,
                    cf.colsubj,
                    cf.colvar_exclude,
                    cf.colvar_atlas,
                    cf.colvar_main,
                    cf.colvar_reg1,
                    cf.colvar_reg2)
    #%% makes file list for subjects
    print('Identifying subjects for atlas creation and hypothesis testing...')
    subTest_fname = []; subTest_IDs = []; subAtlas_fname = []; subAtlas_IDs = []
    for ind in range(len(sub_ID)):
        sub = sub_ID[ind]
        fname = sub_fname[ind]
        if int(subAtlas_idx[ind]) !=1:        
            subTest_fname.append(fname)
            subTest_IDs.append(sub)
        else:
            subAtlas_fname.append(fname)  
            subAtlas_IDs.append(sub)
    len(subTest_IDs)

    count1=0
    subTest_varmain = sp.zeros(len(subTest_IDs)); subTest_varc1 = sp.zeros(len(subTest_IDs)); subTest_varc2 = sp.zeros(len(subTest_IDs))
    for ind in range(len(sub_ID)):
        varmain = reg_var[ind]
        varc1 = reg_cvar1[ind]
        varc2 = reg_cvar2[ind]
        if int(subAtlas_idx[ind]) !=1:  
            subTest_varmain[count1] = varmain
            subTest_varc1[count1] = varc1
            subTest_varc2[count1] = varc2
            count1 +=1

    np.savetxt(cf.out_dir + "/subjects_testing.csv", subTest_IDs, delimiter=",", fmt='%s')
    np.savetxt(cf.out_dir + "/subjects_atlas.csv", subAtlas_IDs, delimiter=",", fmt='%s')
    print(str(len(subAtlas_IDs)) + ' subjects will be used for atlas creation.')
    print(str(len(subTest_IDs)) + ' subjects will be used for hypothesis testing.')
    #%%
    # reads reference data and creates atlas by BrainSync algorithm
    subAtlas_data = load_bfp_data(subAtlas_fname, cf.LenTime)
    if cf.useGroupSync == True:
        print('User Option: Group BrainSync algorithm will be used for atlas creation')
        atlas_data, _, _, _ = groupBrainSync(subAtlas_data)
    else:
        print('User Option: representative subject will be used for atlas creation')
        subRef_data, subRef_num = IDrefsub_BrainSync(subAtlas_data)
        print('Subject number ' + str(subRef_num) + ' will be used for atlas creation')
        atlas_data = generate_avgAtlas(subRef_data, subAtlas_data)
        
    spio.savemat(os.path.join(cf.out_dir + '/atlas.mat'), {'atlas_data': atlas_data})
    del subAtlas_data
    #%% sync and calculates geodesic distances
    subTest_data = load_bfp_data(subTest_fname, cf.LenTime)
    subTest_syndata = sync2atlas(atlas_data, subTest_data)
    subTest_diff = dist2atlas(atlas_data, subTest_syndata)
    spio.savemat(os.path.join(cf.out_dir + '/dist2atlas.mat'), {'subTest_diff': subTest_diff})
    del subTest_data, subTest_syndata
    #%% computes correlation after controlling for two covariates
    rval, pval, pval_fdr = multiLinReg_corr(subTest_diff, subTest_varmain, subTest_varc1, subTest_varc2 )
    #%%
    vis_grayord_sigcorr(pval, rval, cf.outname, cf.out_dir, cf.smooth_iter)
    vis_grayord_sigcorr(pval_fdr, rval, cf.outname + '_fdr', cf.out_dir, cf.smooth_iter)
    visdata_grayord(pval, surf_name, out_dir, cf.smooth_iter, colorbar_lim,colormap)