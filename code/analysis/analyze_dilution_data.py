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
fluct_data = pd.read_csv('../../data/raw_compiled_lineages.csv')
fc_data = pd.read_csv('../../data/raw_compiled_snaps.csv')


# Constants and bounds for size
IP_DIST = 0.065 # In nm / pix
BIRTHFRAC = 0.95 # Fraction of birth lengths to define as "divided" cells.

# Load the stan model
model = mwc.bayes.StanModel('../stan/calibration_factor.stan') #, force_compile=True)

#%%
# Instantiate storage lists for calfactor samples
fc_dfs = []
fluct_dfs = []
for g, d in tqdm.tqdm(fluct_data.groupby(['carbon', 'temp', 'date', 'run_number'])):
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

    # Generate a dataframe with all of the fluctuations
    fluct_df = pd.DataFrame([])
    fluct_df['I_1'] = d['I_1_sub'].values
    fluct_df['I_2'] = d['I_2_sub'].values
    fluct_df['summed'] = d['I_1_sub'].values + d['I_2_sub'].values
    fluct_df['fluct'] = (d['I_1_sub'].values - d['I_2_sub'].values)**2
    fluct_df['alpha_mean'] = samples['alpha'].mean()
    fluct_df['alpha_std'] = np.std(samples['alpha']) 
    fluct_df['carbon'] = g[0]
    fluct_df['temp'] = g[1]
    fluct_df['date'] = g[2]
    fluct_df['run_no'] = g[-1]
    fluct_df['length_1_death'] = d['length_1_death'].values 
    fluct_df['length_1_birth'] = d['length_1_birth'].values 
    fluct_df['length_2_death'] = d['length_2_death'].values 
    fluct_df['length_2_birth'] = d['length_2_birth'].values 
    fluct_df['volume_1_birth'] = d['volume_1_birth'].values
    fluct_df['volume_2_birth'] = d['volume_2_birth'].values
    fluct_df['volume_1_death'] = d['volume_1_death'].values
    fluct_df['volume_2_death'] = d['volume_2_death'].values
    fluct_dfs.append(fluct_df)

    # Compute the fold-change
    _fc = pd.DataFrame([])
    _fc['atc_ngml'] = _fc_data['atc_ngml']
    _fc['date'] = _fc_data['date']
    _fc['run_number'] = _fc_data['run_number']
    _fc['raw_repressors'] = _fc_data['mch_sub'] / samples['alpha'].mean()
    _fc['raw_repressors_max'] =  _fc_data['mch_sub'] / (samples['alpha'].mean() - samples['alpha'].std())
    _fc['raw_repressors_min'] =  _fc_data['mch_sub'] / (samples['alpha'].mean()  + samples['alpha'].std())
    _fc['fold_change'] = _fc_data['yfp_sub'].values / (_fc_data[_fc_data['strain']=='delta']['yfp_sub'].values).mean()
    _fc['yfp_sub'] = _fc_data['yfp_sub']
    _fc['mch_sub'] = _fc_data['mch_sub']
    _fc['area'] = _fc_data['area_pix'] * IP_DIST**2
    _fc['alpha_mu'] = samples['alpha'].mean()
    _fc['alpha_std'] = samples['alpha'].std()
    _fc['carbon'] = g[0]
    _fc['temp'] = g[1]
    _fc['strain'] = _fc_data['strain']
    _fc['length_um'] = _fc_data['length_um']
    _fc['width_um'] = _fc_data['width_um']
    _fc['volume_birth'] = _fc_data['volume_birth']
    _fc['volume_death'] = _fc_data['volume_death']
    fc_dfs.append(_fc)

# Make everythign into one dataframe
fc_df = pd.concat(fc_dfs)
fluct_df = pd.concat(fluct_dfs)

fcs = []
# Make the cell size correction
for g, d in fc_df.groupby(['carbon', 'temp']):
    d = d.copy()
    # Get the birth sizes from the fluctuation data
    growth = fluct_df[(fluct_df['carbon']==g[0]) & (fluct_df['temp']==g[1])]

    # Find the birth length threshold 
    lengths = growth[['length_1_birth', 'length_2_birth']].values.flatten()
    # ind = np.where((np.arange(1, len(lengths)+1)/ len(lengths)) >= BIRTHFRAC)[0][0]
    thresh = np.percentile(lengths, BIRTHFRAC * 100)

     # Designate the cells as "small" or "large"
    d.loc[d['length_um'] < thresh, 'size'] = 'small'
    d.loc[d['length_um'] >= thresh, 'size'] = 'large'

    # Generate the final repressor count
    d['repressors'] = d['raw_repressors']
    d['repressors_max'] = d['raw_repressors_max']
    d['repressors_min'] = d['raw_repressors_min']
    d['threshold'] = thresh
    d.loc[d['size']=='small', 'repressors'] *= 2
    d.loc[d['size']=='small', 'repressors_max'] *= 2
    d.loc[d['size']=='small', 'repressors_min'] *= 2
    d = d[(d['repressors']>0) & (d['fold_change']>=0)]
    fcs.append(d)
fc_df = pd.concat(fcs)

# Concatenate the summary df and save to disk. 
fc_df.to_csv('../../data/analyzed_foldchange.csv', index=False)
fluct_df.to_csv('../../data/analyzed_fluctuations.csv', index=False)

#%%
