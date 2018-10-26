# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import mwc.bayes
import mwc.stats
# TODO: Generate diagnostics document for each run.

# Load in the compiled data
fluct_data = pd.read_csv('../../data/compiled_fluctuations.csv')
fc_data = pd.read_csv('../../data/compiled_fold_change.csv')

# Load the stan model
model = mwc.bayes.StanModel('../stan/hierarchical_dilution_analysis.stan', force_compile=True)

# Instantiate storage lists for calfactor samples
dfs = []
for g, d in fluct_data.groupby(['carbon']):
    d = d.copy()
    print(f'Beginning processing of {g} experiments...')
    # Isolate the fold-change data. 
    _fc_data = fc_data[fc_data['carbon']==g].copy()    
    
    # Compute the mean autofluorescence for each channel. 
    auto = _fc_data[_fc_data['strain']=='auto']
    mean_auto_mch = np.mean(auto['mean_mCherry'] - auto['mCherry_bg_val'])
    mean_auto_yfp = np.mean(auto['mean_yfp'] - auto['yfp_bg_val'])

    # Perform necessary background subtraction for fluctuation measurements. 
    d['I_1_sub'] = (d['I_1'] - d['bg_val'] - mean_auto_mch) * d['area_1']
    d['I_2_sub'] = (d['I_2'] - d['bg_val'] - mean_auto_mch) * d['area_2']

    # Ensure positivity. 
    d = d[(d['I_1_sub'] >= 0) & (d['I_2_sub'] >= 0)]

    # Perform background subtraction for fluorescence measurements. 
    _fc_data['yfp_sub'] = (_fc_data['mean_yfp'] - _fc_data['yfp_bg_val']) * _fc_data['area_pix']
    _fc_data['mch_sub'] = (_fc_data['mean_mCherry'] - _fc_data['mCherry_bg_val']) * _fc_data['area_pix']

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
   
    # Assemble the data dictionary and sample
    data_dict = {'J_exp': np.max(d['rep_idx']),
                 'N_fluct': len(d),
                 'index_1': d['rep_idx'],
                 'J_conc': np.max(mean_fc_data['conc_idx']),
                 'J_conc_exp': np.max(mean_fc_data['rep_idx']),
                 'N_meas': len(mean_fc_data),
                 'index_2': mean_fc_data['conc_idx'],
                 'I_1': d['I_1_sub'],
                 'I_2': d['I_2_sub'],
                 'mcherry': mean_fc_data['mch_sub'],
                 'yfp': mean_fc_data['yfp_sub']} 
    samples = model.sample(data_dict, iter=1000, chains=4, return_df=False,
                           **dict(control=dict(adapt_delta=0.95)))
    
    # Compute the summarized parameters. 
    summary = model.summarize_parameters(parnames=['alpha_1', 'avg_rep', 'fold_change'])
    
    # Convert the dimension to the atc concentration 
    for i, c in enumerate(_fc_data['atc_ngml'].unique()):
        summary.loc[summary['dimension']==int(i+1), 'atc_ngml'] = c
    summary['carbon'] == g 
    dfs.append(summary)
    
    # Save the sampling data frame to disk. 
    samples_df['carbon'] == g
    samples_df.to_csv(f'../../data/{g}_dilution_analysis.csv', index=False)
    
# Concatenate the summary df and save to disk. 
df = pd.concat(dfs)
df.to_csv('../../data/summarized_dilution_analysis.csv', index=False)
    

