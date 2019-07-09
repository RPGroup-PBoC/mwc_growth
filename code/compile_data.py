#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob
import pandas as pd
import sys
import mwc.io

# Define the data directory. 
data_dir = '../code/processing/microscopy/'

# Find all dilution experiments. 
dil_exp = glob.glob(f'{data_dir}*dilution')

fluct_dfs = []
fc_dfs = []
samp_dfs = []
for _, d in enumerate(dil_exp):
    date, run, temp, carbon, op, _ = d.split('/')[-1].split('_')
    info = mwc.io.scrape_frontmatter(f'{d}')
    if info['status'].lower() == 'accepted':
        fluct_df = pd.read_csv(glob.glob(f'{d}/output/{date}_{run}_{temp}_{carbon}_{op}_fluctuations.csv')[0])
        fc_df = pd.read_csv(glob.glob(f'{d}/output/{date}_{run}_{temp}_{carbon}_{op}_foldchange.csv')[0])
        fluct_dfs.append(fluct_df)
        fc_dfs.append(fc_df)

_fluct_df = pd.concat(fluct_dfs, sort=False)
_fluct_df = _fluct_df[['carbon', 'temp', 'date', 'run_no', 'I_1', 'I_2', 'area_1', 'area_2']]
fc_df = pd.concat(fc_dfs, sort=False)
_fluct_df.to_csv('../data/compiled_fluctuations.csv', index=False)
fc_df.to_csv('../data/compiled_fold_change.csv', index=False)
print('all data compiled')
