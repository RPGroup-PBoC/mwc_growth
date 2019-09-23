# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import mwc.bayes
import mwc.stats
import bokeh.io
import bokeh.plotting
import mwc.viz
color, color_list = mwc.viz.bokeh_theme()
import tqdm
bokeh.io.output_notebook()

# Load in the compiled data
fluct_data = pd.read_csv('../../data/compiled_fluctuations.csv')
fc_data = pd.read_csv('../../data/compiled_fold_change.csv')

# Constants and bounds for size
IP_DIST = 0.065 # In nm / pix
min_size = 0.25 / IP_DIST**2 # Number is in µm^2
max_size = 5 / IP_DIST**2 # Number is in µm^2

# Load the stan model
model = mwc.bayes.StanModel('../stan/calibration_factor.stan')

#%%
# Filter the cells on area
fluct_data = fluct_data[(fluct_data['area_1'] >= min_size) & 
                        (fluct_data['area_2'] >= min_size) &
                        (fluct_data['area_1'] <= max_size) & 
                        (fluct_data['area_2'] <= max_size)].copy()

fc_data = fc_data[(fc_data['area_pix'] >= min_size) & 
                 (fc_data['area_pix'] <= max_size)].copy()

# Instantiate storage lists for calfactor samples
fc_dfs = []
fluct_dfs = []
for g, d in tqdm.tqdm(fluct_data.groupby(['carbon', 'temp', 'date', 'run_no'])):
    d = d.copy()
    # Isolate the fold-change data. 
    _fc_data = fc_data[(fc_data['carbon']==g[0]) & (fc_data['temp']==g[1]) & 
    (fc_data['run_number']==g[3]) & (fc_data['date']==g[2])].copy()    
    
    # Compute the mean autofluorescence for each channel. 
    auto = _fc_data[_fc_data['strain']=='auto']
    delta = _fc_data[_fc_data['strain']=='delta']
    mean_auto_mch = np.mean(delta['mean_mCherry'])
    mean_auto_yfp = np.mean(auto['mean_yfp'])

    # Perform necessary background subtraction for fluctuation measurements. 
    d['I_1_sub'] = (d['I_1']- mean_auto_mch) * d['area_1']
    d['I_2_sub'] = (d['I_2']- mean_auto_mch) * d['area_2']

    # Ensure positivity. 
    d = d[(d['I_1_sub'] >= 0) & (d['I_2_sub'] >= 0)]

    # Perform background subtraction for fluorescence measurements. 
    _fc_data['yfp_sub'] = (_fc_data['mean_yfp'] - mean_auto_yfp) * _fc_data['area_pix']
    _fc_data['mch_sub'] = (_fc_data['mean_mCherry'] - mean_auto_mch) * _fc_data['area_pix']

    # Estimate the calibration factor
    fit, samples = model.sample(dict(I1=d['I_1_sub'], I2=d['I_2_sub'], N=len(d)))
    # opt, err = mwc.bayes.estimate_calibration_factor(d['I_1_sub'], d['I_2_sub'])

    # Generate a dataframe with all of the fluctuations
    fluct_df = pd.DataFrame([])
    fluct_df['summed'] = d['I_1_sub'].values + d['I_2_sub'].values
    fluct_df['fluct'] = (d['I_1_sub'].values - d['I_2_sub'].values)**2
    fluct_df['alpha_mean'] = samples['alpha'].mean()
    fluct_df['alpha_std'] = np.std(samples['alpha']) 
    fluct_df['carbon'] = g[0]
    fluct_df['temp'] = g[1]
    fluct_dfs.append(fluct_df)

    # Compute the fold-change
    _fc = pd.DataFrame([])
    _fc['atc_ngml'] = _fc_data['atc_ngml']
    _fc['date'] = _fc_data['date']
    _fc['run_number'] = _fc_data['run_number']
    _fc['repressors'] = _fc_data['mch_sub'] / opt
    _fc['repressors_max'] = _fc_data['mch_sub'] / (opt - err)
    _fc['repressors_min'] = _fc_data['mch_sub'] / (opt + err)
    _fc['fold_change'] = _fc_data['yfp_sub'].values / (_fc_data[_fc_data['strain']=='delta']['yfp_sub'].values).mean()
    _fc['alpha_mu'] = opt
    _fc['alpha_std'] = err
    _fc['carbon'] = g[0]
    _fc['temp'] = g[1]
    _fc['strain'] = _fc_data['strain']
    fc_dfs.append(_fc)


# Concatenate the summary df and save to disk. 
fc_df = pd.concat(fc_dfs)
fc_df.to_csv('../../data/analyzed_foldchange.csv', index=False)
fluct_df = pd.concat(fluct_dfs)
fluct_df.to_csv('../../data/analyzed_fluctuations.csv', index=False)
    

#%%
