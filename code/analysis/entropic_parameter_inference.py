################################################################################
# entropic_parameter_inference.py
# ------------------------------------------------------------------------------
# Author: Griffin Chure
# License: MIT
# 
# Description
# ------------------------------------------------------------------------------
# This script infers the entropic parameters delta S_r and delta S_AI from
# fold-change measurements. It can be run either on all temperature data
# together or individually by toggling the "pooled" variable.
# ##############################################################################
#%%
import numpy as np 
import pandas as pd
import bokeh.io
import mwc.bayes

# Define whether the analysis should be pooled or not. 
pooled = 1
force_compile = True

# Load the fluctuation data
data = pd.read_csv('../../data/analyzed_foldchange.csv', comment='#')

# Keep only the dilution strain and the glucose samples as well as >0 reps
data = mwc.process.condition_filter(data, carbon='glucose')
data['repressors'] *= 1.16
# data = data[data['temp'] != 37]

# Group by date, run_number, and ATC concentration to compute the mean fc
grouped = data.groupby(['date', 'run_number', 'atc_ngml', 'temp']).mean().reset_index()

# Determine the number of unique temperatures and add an identifier. 
grouped['idx'] = grouped.groupby('temp').ngroup() + 1

#%%
# Load the inferential model. 
if pooled == 1:
    model = mwc.bayes.StanModel('../stan/constrained_entropy_estimation.stan', force_compile=force_compile)
else:
    model = mwc.bayes.StanModel('../stan/entropy_estimation.stan', force_compile=force_compile)
#%%
# Assign the data dictionary. 
data_dict = {'ref_temp':37 + 273.15, 'ref_epRA': -13.9, 'ref_epAI':4.5, 'Nns':4.6E6}
if pooled == 1:
    data_dict['J'] =  grouped['idx'].max()
    data_dict['N'] = len(grouped) 
    data_dict['J_REF'] = 2
    data_dict['idx'] = grouped['idx']
    data_dict['temp'] = np.array([32, 37, 42]) + 273.15
    data_dict['repressors'] = grouped['repressors']
    data_dict['foldchange'] = grouped['fold_change']

    # Sample the model.
    fit, samples = model.sample(data_dict, iter=5000, control=dict(adapt_delta=0.99))


    # Extract and summarize the parameters
    params = model.summarize_parameters()

else:
    samps = []
    pars = []
    for g, d in grouped.groupby(['temp']):
        data_dict['N'] = len(d)
        data_dict['temp'] = g + 273.15
        data_dict['repressors'] = d['repressors']
        data_dict['foldchange'] = d['fold_change']
        _, samples = model.sample(data_dict, iter=5000, control=dict(adapt_delta=0.99)) 
        samples['temp'] = g
        samps.append(samples)
        _params = model.summarize_parameters()
        _params['temp'] = g
        pars.append(_params)

    samples = pd.concat(samps)
    params = pd.concat(pars)

#%%
# Rename the dimensions to the correct temperatures. 
keymap = {}
for dim, temp in zip(grouped['idx'].unique(), grouped['temp'].unique()):
    if pooled == 1:
       params.loc[params['dimension']==dim, 'temp'] = temp
    keymap[dim] = temp
    keymap[temp] = dim

# renamed_params = ['epRA_star', 'epAI_star', 'delta_SR', 'delta_SAI', 'sigma']
renamed_params = ['epRA_star', 'delta_SR', 'ref_epRA', 'sigma']
samples_dfs = []
if pooled == 1:
    for g, d in params.groupby(['temp']):
        for p in renamed_params:
            df = pd.DataFrame()
            df['lp__'] = samples['lp__']
            if ((p == 'epRA_star') | (p == 'epAI_star')):
                df['value'] = samples[f'{p}[{keymap[g]}]']   
            else:
                df['value'] = samples[f'{p}']
            df['parameter'] = p
            df['step_id'] = np.arange(len(df)) + 1
            df['temp'] = g
            samples_dfs.append(df)
else:
    for g, d in samples.groupby(['temp']):
        for p in renamed_params:
            df = pd.DataFrame()
            df['value'] = d[f'{p}']
            df['parameter'] = p
            df['carbon'] = 'glucose' 
            df['temp']  = g
            df['step_id'] = np.arange(len(d)) + 1
            samples_dfs.append(df)


# Concatenate the sampling info and save
lf_samples = pd.concat(samples_dfs)
if pooled==1:
    name = 'pooled_entropic_parameter'
else:
    name = 'entropic_parameter'
lf_samples.to_csv(f'../../data/{name}_samples.csv')

# Add the carbon information and save to disk. 
params['carbon'] = 'glucose'
params.to_csv(f'../../data/{name}_summary.csv')


# %%
