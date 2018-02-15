import os
import glob
# Our numerical workhorses
import numpy as np
import pandas as pd
import scipy
# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Seaborn, useful for graphics
import seaborn as sns

# Import the project utils
import sys
sys.path.insert(0, '../../../')

import mwc_growth as mwc

# define variables to use over the script
date = 20170610
username = 'gchure'
run = 'r1'

# list the directory with the data
datadir = '../../../data/mbl_2017/csv/'
files = np.array(os.listdir(datadir))
csv_bool = np.array([str(date) in f and 'csv' in f for f in files])
files = files[np.array(csv_bool)]
# define the patterns in the file names to read them
operator = 'O2'
energy = -13.9
rbs = np.array(['auto', 'delta', 'lacI'])
concentrations = ['000', '002', '003', '004', '006', '008', '010']
# define a dictionary that contains the number of well and the concentration

#===============================================================================
# define the parameter alpha for the automatic gating
alpha = 0.40

# initialize the DataFrame to save the mean expression levels
df = pd.DataFrame()
# read the files and compute the mean YFP value
for i, c in enumerate(concentrations):
    for j, strain in enumerate(rbs):
        # find the file
        try:
            r_file = glob.glob(datadir + str(date) + '_' + \
                    operator + '_' + strain + '_ngATC_' + str(c) + '*csv')
            print(r_file)
            # read the csv file
            dataframe = pd.read_csv(r_file[0])
            # apply an automatic bivariate gaussian gate to the log front
            # and side scattering
            data = mwc.auto_gauss_gate(dataframe, alpha,
                                        x_val='FSC-A', y_val='SSC-A',
                                        log=True)
            # compute the mean and append it to the data frame along the
            # operator and strain
            df = df.append([[date, username, operator, energy,
                        strain, c,
                        data['GFP/YFP-H'].mean(), data['CFP-V1-H'].mean(),
                        data['mCherry-H'].mean()]],
                        ignore_index=True)
        except:
            pass

# rename the columns of the data_frame
df.columns = ['date', 'username', 'operator', 'binding_energy', \
        'rbs', 'ATC', 'mean_YFP_H', 'mean_CFP_H', 'mean_mCherry_H']

# Write to disk
df.to_csv('output/' + str(date) + '_' + run + '_' + operator + \
        '_aTc_titration_MACSQuant.csv', index=False)

#===============================================================================
# Add the comments to the header of the data file
filenames = ['./comments.txt', 'output/' + str(date) + '_' + run + '_' + \
             operator + '_aTc_titration_MACSQuant.csv']
with open('../../../data/' + str(date) + '_' + run + '_' + operator + \
        '_aTc_titration_MACSQuant.csv', 'w') as output:
    for fname in filenames:
        with open(fname) as infile:
            output.write(infile.read())
