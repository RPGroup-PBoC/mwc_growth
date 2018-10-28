# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np 
import pandas as pd 
import tqdm
import mwc.bayes
import mwc.stats

# Load the compiled data
fluct_data = pd.read_csv('../../data/compiled_fluctuations.csv')
fc_data = pd.read_csv('../../data/compiled_fold_change.csv')

# Load the hierarchical model. 
model = mwc.bayes.StanModel('../stan/hierarchical_calibration_factor.stan')

# Loop through each carbon source and analyze independently
_fc_dfs = []
for g, d in fluct_data.groupby('carbon'):
    
    # Isolate the corresponding fold-change data
    _fc_data = fc_data[fc_data['carbon']==g].copy()
    
    # Assign unique identifiers to each replicate. 
    d['index'] = d.groupby(['date', 'run_no']).ngroup() + 1
    
    # Assemble the data dictionary. 
    data_dict = {'J_exp':d['index'].max(),
                'N_fluct': len(d),
                'N_mch': len(_fc_data),
                'index_1': d['index'],
                'I_1': d['I_1'],
                'I_2': d['I_2'],
                'mcherry': _fc_data['total_mCherry']}

    # Sample the model
    samples, samples_df = model.sample(data_dict=data_dict)
    
    # Save the sampling dataframe to disk.
    samples_df.to_csv(f'../../data/{g}_calibration_factor_samples.csv', index=False)
        
    # Extract the mode and hpd of the repressor counts. 
    ind = np.argmax(samples_df['lp__'].values)
    rep_mode = np.zeros(len(_fc_data))
    rep_max = np.zeros(len(_fc_data))
    rep_min = np.zeros(len(_fc_data))
    for i in tqdm.tqdm(range(len(_fc_data)), desc='Counting repressors...'):
        cell = samples_df[f'rep_per_cell[{i+1}]'].values
        rep_mode[i] = cell[ind]
        _min, _max = mwc.stats.compute_hpd(cell, 0.95)
        rep_min[i] = _min
        rep_max[i] = _max
        
    # Add arrays to fold-change data frame
    _fc_data['repressors'] = rep_mode
    _fc_data['repressors_min'] = rep_min
    _fc_data['repressors_max'] = rep_max
    
    # Append the dataframe to the list
    _fc_dfs.append(_fc_data)
fc_df = pd.concat(_fc_dfs)
fc_df.to_csv('../../data/counted_repressors.csv', index=False)
            