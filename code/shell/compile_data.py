#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob
import pandas as pd
import sys
sys.path.insert(0, '../../')
import mwc.io

# Define the data directory. 
data_dir = '../../code/processing/microscopy/'

# Find all dilution experiments. 
dil_exp = glob.glob(f'{data_dir}*dilution')

fluct_dfs = []
fc_dfs = []
for _, d in enumerate(dil_exp):
    # Load the readme file.
    info = mwc.io.scrape_frontmatter(f'{d}')
    
    if info['status'].lower() == 'accepted':
        fluct_df = pd.read_csv(glob.glob(f'{d}/output/*fluctuations.csv')[0])
        fc_df = pd.read_csv(glob.glob(f'{d}/output/*foldchange.csv')[0])
        fluct_dfs.append(fluct_df)
        fc_dfs.append(fc_df)
fluct_df = pd.concat(fluct_dfs)
fc_df = pd.concat(fc_dfs)
fc_df['IPTGuM'] = 0
fluct_df.to_csv('../../data/compiled_fluctuations.csv', index=False)
fc_df.to_csv('../../data/compiled_fold_change.csv', index=False)
print('all data compiled')