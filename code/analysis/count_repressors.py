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
model = mwc.bayes.loadStanModel('../stan/hierarchical_calibration_factor') 

# Instantiate storage lists for calfactor samples
alpha_samples_dfs = []
foldchange_dfs = []
summary_dfs = []

for g, d in fluct_data.groupby(['carbon']):
    d = d.copy()
    print(f'Beginning processing of {g} experiments...')
    # Isolate the fold-change data. 
    _fc_data = fc_data[fc_data['carbon']==g].copy()    
    
    # Compute the mean autofluorescence and delta expression. 
    auto = _fc_data[_fc_data['strain']=='auto']
    delta = _fc_data[_fc_data['strain']=='delta']
    mean_auto_mcherry = np.mean(_fc_data['mean_mCherry'] - _fc_data['mCherry_bg_val'])
    mean_auto_yfp = np.mean(_fc_data['mean_yfp'] - _fc_data['yfp_bg_val'])
    delta_yfp = np.mean((delta['mean_yfp'] - delta['yfp_bg_val']) * delta['area_pix'])
     
    # Subtract background and mean autofluorescence
    d['I_1_sub'] = (d['I_1'] - d['bg_val'] - mean_auto_mcherry) * d['area_1']
    d['I_2_sub'] = (d['I_2'] - d['bg_val'] - mean_auto_mcherry) * d['area_2']
    _fc_data['mCherry_sub'] = (_fc_data['mean_mCherry'] - _fc_data['mCherry_bg_val'] - mean_auto_mcherry) * _fc_data['area_pix']
    
    # Remove negative values
    d = d[(d['I_1_sub'] >= 0) & (d['I_2_sub'] >= 0)].copy()
    _fc_data = _fc_data[_fc_data['mCherry_sub'] >= 0].copy()
     
    # Add identifiers. 
    d['idx'] = d.groupby(['date', 'run_no']).ngroup() + 1
     
    # Assemble the data dictionary
    data_dict = {'J_exp':d['idx'].max(),
                'N_fluct': len(d),
                'N_mch': len(_fc_data),
                'index_1':d['idx'],
                'I_1': d['I_1_sub'],
                'I_2': d['I_2_sub'],
                'mcherry': _fc_data['mCherry_sub']}
    
    # Sample the model. 
    print('sampling...')
    samples = model.sampling(data_dict)
    print('finished sampling!')
    
    # Convert to a dataframe and extract the important parameters. 
    samples_df = samples.to_dataframe()
    samples_df['carbon'] = g
    alpha_samples_dfs.append(samples_df)
    
    # Extract the repressor counts
    fit = samples.extract()
    min_alpha, max_alpha = mwc.stats.compute_hpd(fit['alpha_1'], 0.95)
    med_alpha = np.median(fit['alpha_1'])
    min_rep = np.empty(len(_fc_data))
    max_rep = np.empty(len(_fc_data))
    med_rep = np.empty(len(_fc_data)) 
    for i in range(len(_fc_data)): 
        hpd_min, hpd_max = mwc.stats.compute_hpd(fit['rep_per_cell'][i], 0.95)
        _med_rep = np.median(fit['rep_per_cell'][i])
        min_rep[i] = hpd_min
        max_rep[i] = hpd_max
        med_rep[i] = _med_rep
   
    # Insert the repressor counts into the dataframe. 
    _fc_data['min_rep'] = min_rep
    _fc_data['max_rep'] = max_rep
    _fc_data['median_rep'] = med_rep
    _fc_data['min_alpha'] = min_alpha
    _fc_data['max_alpha'] = max_alpha
    _fc_data['median_alpha'] = med_alpha
    
    # Recompute the fold-change. 
    _fc_data['fold_change'] = (_fc_data['mean_yfp'] - _fc_data['yfp_bg_val'] - mean_auto_yfp) * _fc_data['area_pix'] /delta_yfp 
   
    # Add the longform foldchange data to the storage list
    foldchange_dfs.append(_fc_data)
    
    # Compute the summary statistics. 
    _grouped = _fc_data.groupby(['strain', 'atc_ngml'])[['min_rep', 'max_rep', 'median_rep', 'fold_change']].mean().reset_index()
    summary_dfs.append(_grouped)

# Assemble the total dataframes and compute the summary. 
alpha_samples_df = pd.concat(alpha_samples_dfs)
foldchange_df = pd.concat(foldchange_dfs)
summary_df = pd.concat(summary_dfs)

# Save all to disk. 
print('Saving data to disk...')
alpha_samples_df.to_csv('../../data/calibration_factor_samples.csv', index=False)
foldchange_df.to_csv('../../data/foldchange_repressors.csv', index=False)
summary_df.to_csv('../../data/foldchange_summary.csv', index=False)
print('Completed. Thank you come again.')
    
