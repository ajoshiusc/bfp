import csv
import os
from tqdm import tqdm
import scipy as sp
import numpy as np
import scipy.io as spio


def read_oasis3_data(csv_fname,
                       data_dir,
                       reg_var_name='Age',
                       num_sub=5,
                       reg_var_positive=1,
                       len_time=20,
                       data_field='dtseries'):
    """ reads fcon1000 csv and data"""

    count1 = 0
    sub_ids = []
    reg_var = []
    pbar = tqdm(total=num_sub)

    with open(csv_fname, newline='') as csvfile:
        creader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
        for row in creader:

            # read the regression variable
            rvar = row[reg_var_name]

            # Read the filtered data by default
            fname = row['FileName']

            # If the data does not exist for this subject then skip it
            if not os.path.isfile(fname):
                continue
            num_v_t = spio.loadmat(fname)[data_field].shape
    
            # Check if there are enough time points
            if num_v_t[1]<len_time:
                print(fname + ' doesn\'t have enough timepoints' + str(num_v_t[1]) + '/' + str(len_time))
                continue

            if len(rvar)==0 or (reg_var_positive == 1 and np.float64(rvar) < 0):
                continue

            if count1 == 0:
                sub_data_files = []

            # Truncate the data at a given number of time samples This is needed because
            # BrainSync needs same number of time sampples
            sub_data_files.append(fname)
            sub_ids.append(row['Subject'])
            reg_var.append(float(rvar))

            count1 += 1
            pbar.update(1)  # update the progress bar
            #print('%d,' % count1, end='')
            if count1 == num_sub:
                break

    pbar.close()
    print('CSV file and the data has been read\nThere are %d subjects' %
          (len(sub_ids)))

    return sub_ids, sp.array(reg_var), sub_data_files
