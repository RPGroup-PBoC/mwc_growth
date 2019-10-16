#%%
import numpy as np 
import pandas as pd
import bokeh.io
import bokeh.plotting
import mwc.viz
import mwc.bayes
colors, color_list = mwc.viz.bokeh_theme()
bokeh.io.output_notebook()

# Load the fluctuation data
_data = pd.read_csv('../../data/analyzed_foldchange.csv')

# Keep only the dilution strain and the glucose samples as well as >0 reps
data = _data[(_data['strain']=='dilution') &
            (_data['carbon']=='glucose') &
            (_data['repressors'] > 10)].copy()

# Group by date, run_number, and ATC concentration to compute the mean fc
grouped = data.groupby(['date', 'run_number', 'atc_ngml']).mean().reset_index()

# Determine the number of unique temperatures and add an identifier. 
grouped['idx'] = grouped.groupby('temp').ngroup() + 1

# Load the inferential model. 
model = mwc.bayes.StanModel('../stan/entropy_estimation.stan', force_compile=True)
#%%
# Assign the data dictionary. 
data_dict = {'J': grouped['idx'].max(), 
             'N':len(grouped), 
             'idx': grouped['idx'],
             'ref_temp': 37,
             'ref_epRA': -13.9,
             'ref_epAI': 4.5,
             'Nns': 4.6E6,
             'temp': grouped['temp'].values,
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

renamed_params = ['delta_S_DNA', 'sigma', 'true_epRA']
samples_dfs = []
for k, v in keymap.items():
    for p in renamed_params:
        df = pd.DataFrame()
        if (p == 'delta_S_DNA') | (p=='delta_S_ALLO'):
            df['value'] = samples[f'{p}']
        elif (p == 'true_epRA'):
            df['value'] = samples[f'{p}']
        else:
            df['value'] = samples[f'{p}[{k}]']
        df['lp__'] = samples['lp__']
        df['parameter'] = p
        df['temp'] = v
        df['carbon'] = 'glucose'
        samples_dfs.append(df)

# Concatenate the sampling info and save
lf_samples = pd.concat(samples_dfs)
lf_samples.to_csv('../../data/entropic_parameter_samples.csv')

# Add the carbon information and save to dis,. 
params['carbon'] = 'glucose'
params.to_csv('../../data/entropic_parameter_summary.csv')


#%%


#%%
