""" This module contains helpful utility function for running statistics using BFP """

import csv
import os
import scipy as sp
import scipy.io as spio
from brainsync import normalizeData, brainSync
from tqdm import tqdm


def read_fcon1000_data(csv_fname,
                       data_dir,
                       reg_var_name='Verbal IQ',
                       num_sub=5,
                       len_time=250):
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
            fname = os.path.join(
                data_dir, row['ScanDir ID'] + '_rest_bold.32k.GOrd.filt.mat')

            # If the data does not exist for this subject then skip it
            if not os.path.isfile(fname) or int(row['QC_Rest_1']) != 1:
                continue

            # Load data and normalize it
            data = spio.loadmat(fname)
            data = data['dtseries'].T
            data, _, _ = normalizeData(data)

            if count1 == 0:
                sub_data = sp.zeros((len_time, data.shape[1], num_sub))

            # Truncate the data at a given number of time samples This is needed because
            # BrainSync needs same number of time sampples
            sub_data[:, :, count1] = data[:len_time, ]
            sub_ids.append(row['ScanDir ID'])
            reg_var.append(float(rvar))

            count1 += 1
            pbar.update(1) # update the progress bar
            #print('%d,' % count1, end='')
            if count1 == num_sub:
                break

    pbar.close()
    print('CSV file and the data has been read\nThere are %d subjects' %
          (len(sub_ids)))

    return sub_ids, reg_var, sub_data


def sync2atlas(atlas, sub_data):
    print('Syncing to atlas, assume that the data is normalized')

    for ind in tqdm(range(sub_data.shape[2])):
        sub_data[:, :, ind], _ = brainSync(X=atlas, Y=sub_data[:, :, ind])
