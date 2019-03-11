#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob
import pandas as pd
import sys
sys.path.insert(0, '../')
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
        fluct_df = pd.read_csv(glob.glob(f'{d}/output/*fluctuations.csv')[0])
        fc_df = pd.read_csv(glob.glob(f'{d}/output/*foldchange.csv')[0])
        # Remove a problem data point.
        if (date == '20181026') & (carbon == 'glycerol'):
            fc_df = fc_df[fc_df['atc_ngml'] != 0.7] 
        fluct_dfs.append(fluct_df)
        fc_dfs.append(fc_df)
_fluct_df = pd.concat(fluct_dfs, sort=False)
fc_df = pd.concat(fc_dfs, sort=False)

# Compute the background subtraction using the mean autofluorescence. 
fluct_dfs = []
fc_dfs = []
for g, d in _fluct_df.groupby(['carbon', 'date', 'run_no']):
    d = d.copy()
    _fc_data = fc_df[(fc_df['carbon']==g[0]) &
                     (fc_df['date']==g[1]) & 
                     (fc_df['run_number']==g[2])].copy()
#     auto = _fc_data[_fc_data['strain']=='auto']
#     mean_auto_mch = auto['mean_mCherry'].mean()
#     mean_auto_yfp = auto['mean_yfp'].mean()
    
#     # Subtract from the fold-change data
#     _fc_data['mean_mCherry'] -= mean_auto_mch
#     _fc_data['mean_yfp'] -= mean_auto_yfp 
#     _fc_data['total_mCherry'] = _fc_data['mean_mCherry'] * _fc_data['area_pix']
    
#     # Compute the repressors per cell from day-to-day calibration factor. 
#     _fc_data['repressors'] = _fc_data['total_mCherry'] / _fc_data['alpha_mean']
#     _fc_data['repressors_min'] = _fc_data['total_mCherry'] / _fc_data['alpha_hpd_min']
#     _fc_data['repressors_max'] = _fc_data['total_mCherry'] / _fc_data['alpha_hpd_max']
     
#     # Subtract from the fluctuation data. 
#     d['I_1'] = (d['I_1'] - mean_auto_mch) * d['area_1']
#     d['I_2'] = (d['I_2'] - mean_auto_mch) * d['area_2']
    
#     # Ensure positivity
#     d = d[(d['I_1'] >= 0) & (d['I_2'] >= 0)]
    
#     # Compute the fluctuations and squared differences for simple plotting. 
#     d['summed'] = d['I_1'] + d['I_2']
#     d['sq_fluct'] = (d['I_1'] - d['I_2'])**2

#     # Append the dataframes to lists and concatenate.
    fluct_dfs.append(d)
    fc_dfs.append(_fc_data)
    
fluct_df = pd.concat(fluct_dfs, sort=False)
fc_df = pd.concat(fc_dfs, sort=False)
fluct_df.to_csv('../data/compiled_fluctuations.csv', index=False)
fc_df.to_csv('../data/compiled_fold_change.csv', index=False)


# Find all microscopy growth experiments. 
data_dir = '../code/processing/growth_curves_microscopy/'
growth_exp = glob.glob(f'{data_dir}*growth')
growth_dfs = []
for _, d in enumerate(growth_exp):
    # Load the readme file.
    info = mwc.io.scrape_frontmatter(f'{d}')
    
    if info['status'].lower() == 'accepted':
        growth_df = pd.read_csv(glob.glob(f'{d}/output/*growth.csv')[0])
        growth_dfs.append(growth_df)
growth_df = pd.concat(growth_dfs)
growth_df.to_csv('../data/compiled_growth_microscopy.csv', index=False)
print('all data compiled')