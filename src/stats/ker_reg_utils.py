
from brainsync import normalizeData
from sklearn.kernel_ridge import KernelRidge as KRR
import itertools
import numpy as np
from stats_utils import pair_dist, pair_dist_simulation
import scipy.io as spio
from tqdm import tqdm
from functools import partial
from multiprocessing import Pool



def kernel_regression_choose_gamma(bfp_path,
                      sub_files,
                      reg_var,
                      nperm=1000,
                      len_time=235,
                      num_proc=4,
                      fdr_test=False):
    """ Choose gamma for Kernel Regression """



    # Normalize the variable
    reg_var, _, _ = normalizeData(reg_var)

    # Get the number of vertices from a file
    num_vert = spio.loadmat(sub_files[0])['dtseries'].shape[0]
    num_sub = len(sub_files)
    pairs = np.array(list(itertools.combinations(range(num_sub), r=2)))
    num_pairs = len(pairs)

    fmri_diff = np.zeros((num_vert, num_pairs))
    regvar_diff = np.zeros(num_pairs)
    # added for simulation
    labs = spio.loadmat(
        '/ImagePTE1/ajoshi/code_farm/bfp/supp_data/USCLobes_grayordinate_labels.mat')['labels']
    roi = (labs == 200)  # R. Parietal Lobe

    pairdistfunc = pair_dist_simulation

    if num_proc > 1:
        pool = Pool(num_proc)

        results = pool.imap(
            partial(pairdistfunc,
                    sub_files=sub_files,
                    reg_var=reg_var,
                    len_time=len_time,
                    roi=roi), pairs)

        ind = 0
        for res in results:
            fmri_diff[:, ind] = res[0]
            regvar_diff[ind] = res[1]
            ind += 1

    else:
        for ind in tqdm(range(len(pairs))):

            fmri_diff[:, ind], regvar_diff[ind] = pairdistfunc(
                sub_files=sub_files,
                reg_var=reg_var,
                len_time=len_time,
                rand_pair=pairs[ind],roi=roi)

    kr = KRR(kernel='precomputed') #, alpha=1.1)
    D = np.zeros((num_sub, num_sub))
    pval_kr = np.zeros(num_vert)
    #5  # checked by brute force #5 gives a lot of significance  # bandwidth for RBF

    nperm = 50

    rho = np.zeros(num_vert)
    num_sub_val = 5
    gamma_values = np.arange(1e-8,15,.1)
    rho_all=np.zeros(len(gamma_values))

    roi_ind, _ = np.where(roi)

    for i, gamma in enumerate(gamma_values):
        for v in roi_ind[::5]: #range(0,num_vert,100):
            D = np.zeros((num_sub, num_sub))
            D[pairs[:, 0], pairs[:, 1]] = fmri_diff[v, :]

            D = D+D.T  # make it symmetric

            D = np.exp(-gamma * D)
            # Do this in a split train test split

            D_train = D[:num_sub-num_sub_val,:num_sub-num_sub_val]
            kr.fit(D_train, reg_var[:num_sub-num_sub_val])

            D_val = D[num_sub-num_sub_val:,:num_sub-num_sub_val]
            pred_v = kr.predict(D_val)

            if np.var(pred_v)<1e-6:
                rho[v] = 0
            else:
                rho[v] = np.corrcoef(pred_v,reg_var[num_sub-num_sub_val:])[0,1]

        rho_all[i]=np.mean(rho)
        

        print(gamma,np.mean(rho))

    print(np.argmax(rho_all), gamma_values[np.argmax(rho_all)])

    return gamma_values[np.argmax(rho_all)]


