import csv
import os
from tqdm import tqdm



def read_oasis3_data(csv_fname,
                       data_dir,
                       reg_var_name='Verbal IQ',
                       num_sub=5,
                       reg_var_positive=1,
                       gord=1):
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

            if gord == 1:
                # Read the filtered data by default
                subname = row['Subject']
                fname = os.path.join(data_dir, subname,'func',subname+'_rest_bold.32k.GOrd.filt.mat')
                #fname = os.path.join(
                #    data_dir, row['ScanDir ID'] + '_rest_bold.32k.GOrd.mat')
            else:
                fname = os.path.join(
                    data_dir, row['ScanDir ID'] + '_rest_bold.BOrd.mat')

            # If the data does not exist for this subject then skip it
            if not os.path.isfile(fname):
                continue

            if reg_var_positive == 1 and sp.float64(rvar) < 0:
                continue

            if count1 == 0:
                sub_data_files = []

            # Truncate the data at a given number of time samples This is needed because
            # BrainSync needs same number of time sampples
            sub_data_files.append(fname)
            sub_ids.append(row['ScanDir ID'])
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
