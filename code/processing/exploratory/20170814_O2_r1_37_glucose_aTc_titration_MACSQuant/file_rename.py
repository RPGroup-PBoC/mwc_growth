"""
filename: file_rename.py
date: 2017-08-25
author: Griffin Chure
purpose: Rename FCS files to contain strain and ATC concentration
information
"""
import numpy as np
import pandas as pd
import glob
import os
import shutil

# Define the parameters for the specific experiment.
DATE = '2017-08-14'
OPERATOR = 'O2'
TEMP = 37
CARBON = 'glucose'
RUN_NO = 'r1'

# ----------------------------------------------------------------------------
# Shouldn't need to adjust anything below here
# ----------------------------------------------------------------------------

# Set up the dictionary for well numbers.
atc_concs = [0, 0, 2, 0, 4, 6, 8, 10]
IPTG_concs = [0, 0.1, 5, 10, 25, 50, 75, 100, 250, 500, 1000, 5000]
strain_list = ['auto', 'delta']
strain_list += ['dilution'] * (8 - len(strain_list))
ATC_WELLS = {}
STRAIN_WELLS = {}
IPTG_WELLS = {}
window = np.arange(1, 96, 8)
for i, row in enumerate(window):
    for j, well in enumerate(range(row, row+8)):
        if well < 10:
            prefix = '000'
        else:
            prefix = '00'
        ATC_WELLS[prefix + str(well)] = atc_concs[j]
        STRAIN_WELLS[prefix + str(well)] = strain_list[j]
        IPTG_WELLS[prefix + str(well)] = IPTG_concs[i]

# See if the destination directory exists. If not, then make it.
target = '../../../data/flow/fcs/{0}/renamed/'.format(DATE)
if os.path.isdir(target) is False:
    print('making new directory.')
    os.mkdir(target)


# Grab all of the files with the matching date.
files = glob.glob('../../../data/flow/fcs/{0}/*.fcs'.format(DATE))
for f in files:
    well = f.split('.')[-2]
    atc = ATC_WELLS[well]
    iptg = IPTG_WELLS[well]
    strain = STRAIN_WELLS[well]
    # Define the new file name.
    _DATE = "".join(DATE.split('-'))
    new_name = "{0}_{1}_{2}_{3}_{4}_{5}_{6}ngml_aTc_{7}uM_IPTG.fcs".format(_DATE, OPERATOR, TEMP, CARBON, RUN_NO, strain, atc, iptg)
    shutil.copy(f, target + new_name)
    print(f.split('/')[-1] + ' --> ' + new_name)
