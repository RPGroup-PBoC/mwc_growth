# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import mwc.bayes
import mwc.stats
import tqdm

# Load in the compiled data
fluct_data = pd.read_csv('../../data/compiled_fluctuations.csv')
fc_data = pd.read_csv('../../data/compiled_fold_change.csv')

# Instantiate storage lists for calfactor samples
fc_dfs = []
fluct_dfs = []
for g, d in tqdm.tqdm(fluct_data.groupby(['carbon', 'temp', 'date', 'run_no'])):
    d = d.copy()
    # Isolate the fold-change data. 
    _fc_data = fc_data[(fc_data['carbon']==g[0]) & (fc_data['temp']==g[1]) &
                      (fc_data['date']==g[2]) & (fc_data['run_number']==g[-1])].copy()    
    
    # Compute the mean autofluorescence for each channel. 
    auto = _fc_data[_fc_data['strain']=='auto']
    delta = _fc_data[_fc_data['strain']=='delta']
    mean_auto_mch = np.median(auto['mean_mCherry'])
    mean_auto_yfp = np.median(auto['mean_yfp'])

    # Perform necessary background subtraction for fluctuation measurements. 
    d['I_1_sub'] = (d['I_1']- mean_auto_mch) * d['area_1']
    d['I_2_sub'] = (d['I_2']- mean_auto_mch) * d['area_2']

    # Ensure positivity. 
    d = d[(d['I_1_sub'] >= 0) & (d['I_2_sub'] >= 0)]

    # Perform background subtraction for fluorescence measurements. 
    _fc_data['yfp_sub'] = (_fc_data['mean_yfp'] - mean_auto_yfp) * _fc_data['area_pix']
    _fc_data['mch_sub'] = (_fc_data['mean_mCherry'] - mean_auto_mch) * _fc_data['area_pix']

    # Add relevant identifiers to each data set. 
    d['date_idx'] = d.groupby(['date']).ngroup() + 1
    d['rep_idx'] = d.groupby(['date', 'run_no']).ngroup() + 1

    # Compute the means for the fold-change data
    mean_fc_data = _fc_data.groupby(['date', 'atc_ngml', 'run_number', 'strain']).mean().reset_index()

    # Enforce positivity of YFP and mCherry.
    mean_fc_data = mean_fc_data[mean_fc_data['strain'] != 'auto'].copy()
    mean_fc_data = mean_fc_data[(mean_fc_data['yfp_sub'] >= 0) & (mean_fc_data['mch_sub'] >= 0)]

    # Add appropriate identifiers to mean data.
    mean_fc_data['conc_idx'] = mean_fc_data.groupby(['atc_ngml', 'strain']).ngroup() + 1
    mean_fc_data['rep_idx'] = mean_fc_data.groupby(['date', 'atc_ngml', 'run_number', 'strain']).ngroup() + 1
   
    # Estimate teh calibration factor
    opt, err = mwc.bayes.estimate_calibration_factor(d['I_1_sub'], d['I_2_sub'])

    # Generate a dataframe with all of the fluctuations
    fluct_df = pd.DataFrame([])
    fluct_df['summed'] = d['I_1_sub'].values + d['I_2_sub'].values
    fluct_df['fluct'] = (d['I_1_sub'].values - d['I_2_sub'].values)**2
    fluct_df['alpha_mu'] = opt
    fluct_df['alpha_std'] = err
    fluct_df['carbon'] = g[0]
    fluct_df['temp'] = g[1]
    fluct_dfs.append(fluct_df)

    # Compute the fold-change
    _fc = pd.DataFrame([])
    _fc['repressors'] = _fc_data['mch_sub'] / opt
    _fc['repressors_max'] = _fc_data['mch_sub'] / (opt - err)
    _fc['repressors_min'] = _fc_data['mch_sub'] / (opt + err)
    _fc['fold_change'] = _fc_data['yfp_sub'].values / (_fc_data[_fc_data['strain']=='delta']['yfp_sub'].values).mean()
    _fc['alpha_mu'] = opt
    _fc['alpha_std'] = err
    _fc['carbon'] = g[0]
    _fc['temp'] = g[1]
    fc_dfs.append(_fc)
    
# Concatenate the summary df and save to disk. 
fc_df = pd.concat(fc_dfs)
fc_df.to_csv('../../data/analyzed_foldchange.csv', index=False)
fluct_df = pd.concat(fluct_dfs)
fluct_df.to_csv('../../data/analyzed_fluctuations.csv', index=False)
    


#%%
