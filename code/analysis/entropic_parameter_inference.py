#%%
import numpy as np 
import pandas as pd
import bokeh.io
import mwc.bayes

# Define whether the analysis should be pooled or not. 
pooled = 1
force_compile = True

# Load the fluctuation data
data = pd.read_csv('../../data/analyzed_foldchange.csv')

# Keep only the dilution strain and the glucose samples as well as >0 reps
data = data[(data['strain']=='dilution') &
            (data['carbon']=='glucose') & (data['fold_change'] >= 0) &
            (data['repressors'] > 0) & (data['temp'] != 37)].copy()

# Group by date, run_number, and ATC concentration to compute the mean fc
grouped = data.groupby(['date', 'run_number', 'atc_ngml', 'temp']).mean().reset_index()

# Determine the number of unique temperatures and add an identifier. 
grouped['idx'] = grouped.groupby('temp').ngroup() + 1

#%%
# Load the inferential model. 
if pooled == 1:
    model = mwc.bayes.StanModel('../stan/pooled_entropy_estimation.stan', force_compile=force_compile)

else:
    model = mwc.bayes.StanModel('../stan/entropy_estimation.stan', force_compile=force_compile)
#%%
# Assign the data dictionary. 
data_dict = {'J': grouped['idx'].max(), 
             'N':len(grouped), 
             'idx': grouped['idx'],
             'ref_temp': 37 + 273.15,
             'ref_epRA': -13.9,
             'ref_epAI': 4.5,
             'Nns': 4.6E6,
             'temp': grouped['temp'].unique() + 273.15,
             'repressors': grouped['repressors'],
             'foldchange': grouped['fold_change']}

# Sample the model.
fit, samples = model.sample(data_dict, iter=5000, control=dict(adapt_delta=0.99))

#%%
# Extract and summarize the parameters
params = model.summarize_parameters()
params

#%%
# Rename the dimensions to the correct temperatures. 
keymap = {}

for dim, temp in zip(grouped['idx'].unique(), grouped['temp'].unique()):
    params.loc[params['dimension']==dim, 'temp'] = temp
    keymap[dim] = temp

renamed_params = ['true_epRA', 'true_epAI', 'epRA_star', 'epAI_star', 'delta_S', 'delta_S_vib', 'sigma']
samples_dfs = []
for k, v in keymap.items():
    for p in renamed_params:
        df = pd.DataFrame()
        if pooled == 1:
            if ('true' in p) | (p == 'delta_S') | (p == 'sigma') | (p == 'delta_S_vib'):
                df['value'] = samples[f'{p}']
            else:
                df['value'] = samples[f'{p}[{k}]']
        else:
            df['value'] = samples[f'{p}[{k}]']
        df['lp__'] = samples['lp__']
        df['parameter'] = p
        df['temp'] = v
        df['carbon'] = 'glucose'
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
params.to_csv('../../data/{name}_summary.csv')




#%%


#%%
