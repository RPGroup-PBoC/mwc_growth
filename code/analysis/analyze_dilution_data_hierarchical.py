# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import mwc.bayes
import mwc.stats
import joblib
import tqdm

# Load in the compiled data
fluct_data = pd.read_csv('../../data/compiled_fluctuations.csv')
fc_data = pd.read_csv('../../data/compiled_fold_change.csv')

# Apply an area filter
# Constants and bounds for size
IP_DIST = 0.065 # In nm / pix
min_size = 0.5 / IP_DIST**2 # Number is in µm 
max_size = 4 / IP_DIST**2 # Number is in µm

fluct_data = fluct_data[(fluct_data['area_1'] >= min_size) & 
                        (fluct_data['area_2'] >= min_size) &
                        (fluct_data['area_1'] <= max_size) & 
                        (fluct_data['area_2'] <= max_size)].copy()

fc_data = fc_data[(fc_data['area_pix'] >= min_size) & 
                 (fc_data['area_pix'] <= max_size)]

# Instantiate storage lists for calfactor samples
fc_dfs = []
fluct_dfs = []

# Load the inferential model
model = mwc.bayes.StanModel('../stan/hierarchical_calibration_factor.stan')

#%%
# Perform inference for each unique sample
for g, d in tqdm.tqdm(fluct_data.groupby(['carbon', 'temp'])):
    d = d.copy()
    d['idx'] = d.groupby(['date']).ngroup() + 1
    # d['idx'] = 1

    # Isolate the fold-change data. 
    _fc_data = fc_data[(fc_data['carbon']==g[0]) & (fc_data['temp']==g[1])].copy()    
    _fc_data['idx'] = _fc_data.groupby(['date', 'run_number']).ngroup() + 1
    
    # Compute the mean autofluorescence for each channel. 
    auto = _fc_data[_fc_data['strain']=='auto']
    delta = _fc_data[_fc_data['strain']=='delta']
    mean_auto_mch = np.median(auto['mean_mCherry'])
    mean_auto_yfp = np.median(auto['mean_yfp'])

    # Match up the autofluorescence values for the particular days. 
    auto_mch_dict = {}
    auto_yfp_dict = {}
    for _g, _ in d.groupby(['date', 'run_no', 'idx']):  
        delta = _fc_data[(_fc_data['date']==_g[0]) & (_fc_data['run_number']==_g[1]) 
                        & (_fc_data['strain']=='delta')]
        auto = _fc_data[(_fc_data['date']==_g[0]) & (_fc_data['run_number']==_g[1]) 
                        & (_fc_data['strain']=='auto')]
        _fc_data.loc[(_fc_data['date']==_g[0]) & (_fc_data['run_number']==_g[1]), 'idx'] = _g[-1]
        _fc_data.loc[(_fc_data['date']==_g[0]) &\
                 (_fc_data['run_number']==_g[1]),
                 'mean_auto_mch'] = delta['mean_mCherry'].mean()
        _fc_data.loc[(_fc_data['date']==_g[0]) &\
                 (_fc_data['run_number']==_g[1]),
                 'mean_auto_yfp'] = auto['mean_yfp'].mean()
        d.loc[(d['date']==_g[0]) & (d['run_no']==_g[1]), 'mean_auto'] = mean_auto_mch

    # Perform necessary background subtraction for fluctuation measurements. 
    d['I_1_sub'] = (d['I_1']- mean_auto_mch) * d['area_1']
    d['I_2_sub'] = (d['I_2']- mean_auto_mch) * d['area_2']

    # Ensure positivity. 
    d = d[(d['I_1_sub'] >= 0) & (d['I_2_sub'] >= 0)]

    # Assemble the data dictionary for sampling
    data_dict = dict(J_exp=d['idx'].max(), N_fluct=len(d), index_1=d['idx'],
                    I_1=d['I_1_sub'], I_2=d['I_2_sub'])
    _, samples = model.sample(data_dict, iter=2000)

    # Compute the statistics from the sampling output
    parnames = [f'alpha_2[{i}]' for i in d['idx'].unique()]
    parnames.append('alpha_1')
    params = mwc.stats.compute_statistics(samples, varnames=parnames, logprob_name='lp__')


    # Generate a mapping of number to calibration factor. 
    median_idx = {i:params[(params['parameter']==f'alpha_2[{i}]')]['median'] for i in d['idx'].unique()} 
    mean_idx = {i:params[(params['parameter']==f'alpha_2[{i}]')]['mean'] for i in d['idx'].unique()} 
    mode_idx = {i:params[(params['parameter']==f'alpha_2[{i}]')]['mode'] for i in d['idx'].unique()} 
    min_idx = {i:params[(params['parameter']==f'alpha_2[{i}]')]['hpd_max'] for i in d['idx'].unique()} 
    max_idx = {i:params[(params['parameter']==f'alpha_2[{i}]')]['hpd_min'] for i in d['idx'].unique()}

    # Generate a dataframe with all of the fluctuations
    fluct_df = pd.DataFrame([])
    fluct_df['summed'] = d['I_1_sub'].values + d['I_2_sub'].values
    fluct_df['fluct'] = (d['I_1_sub'].values - d['I_2_sub'].values)**2
    fluct_df['idx'] = d['idx'].values
    fluct_df['hyper_alpha_median'] = params[params['parameter']=='alpha_1']['median']
    fluct_df['hyper_alpha_mean'] = params[params['parameter']=='alpha_1']['mean']
    fluct_df['hyper_alpha_mode'] = params[params['parameter']=='alpha_1']['mode']
    fluct_df['hyper_alpha_min'] = params[params['parameter']=='alpha_1']['hpd_min']
    fluct_df['hyper_alpha_max'] = params[params['parameter']=='alpha_1']['hpd_max']
    fluct_df['carbon'] = g[0]
    fluct_df['temp'] = g[1]
    for _g, _ in fluct_df.groupby(['idx']):
        fluct_df.loc[fluct_df['idx']==_g, 'alpha_median'] = median_idx[_g] 
        fluct_df.loc[fluct_df['idx']==_g, 'alpha_mu'] = mean_idx[_g] 
        fluct_df.loc[fluct_df['idx']==_g, 'alpha_mode'] = mode_idx[_g] 
        fluct_df.loc[fluct_df['idx']==_g, 'alpha_min'] = min_idx[_g] 
        fluct_df.loc[fluct_df['idx']==_g, 'alpha_max'] = max_idx[_g] 

    # Perform background subtraction for fluorescence measurements. 
    _fc_data['yfp_sub'] = (_fc_data['mean_yfp'] - _fc_data['mean_auto_yfp']) * _fc_data['area_pix']
    _fc_data['mch_sub'] = (_fc_data['mean_mCherry'] - _fc_data['mean_auto_mch']) * _fc_data['area_pix']

    # Compute the means for the fold-change data
    mean_fc_data = _fc_data.groupby(['date', 'atc_ngml', 'run_number', 'strain']).mean().reset_index()

    # Enforce positivity of YFP and mCherry.
    mean_fc_data = mean_fc_data[mean_fc_data['strain'] != 'auto'].copy()
    mean_fc_data = mean_fc_data[(mean_fc_data['yfp_sub'] >= 0) & (mean_fc_data['mch_sub'] >= 0)]
   
    # Compute the fold-change
    _fc = _fc_data.copy()
    _fc['repressors'] = _fc_data['mch_sub'] / median_idx[_fc['idx'].values]
    _fc['repressors_max'] = _fc_data['mch_sub'] / min_idx[_fc['idx'].values]
    _fc['repressors_min'] = _fc_data['mch_sub'] / max_idx[_fc['idx'].values]
    for _g, _ in _fc.groupby(['idx']):
        _delta = d[d['strain']=='delta']['yfp_sub'].mean()
        _fc.loc[_fc['idx']==_g, 'fold_change'] = _fc_data['yfp_sub'].values / _delta
    
    _fc['carbon'] = g[0]
    _fc['temp'] = g[1]
    fc_dfs.append(_fc) 
    fluct_dfs.append(fluct_df)

# Concatenate the summary df and save to disk. 
fc_df = pd.concat(fc_dfs)
fc_df.to_csv('../../data/analyzed_foldchange_hierarchical.csv', index=False)
fluct_df = pd.concat(fluct_dfs)
fluct_df.to_csv('../../data/analyzed_fluctuations_hierarchical.csv', index=False)
#%%
