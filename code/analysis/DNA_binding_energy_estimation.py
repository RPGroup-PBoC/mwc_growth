# %%
import numpy as np
import pandas as pd
import mwc.bayes
import tqdm

#  Load the data.
data = pd.read_csv('../../data/analyzed_foldchange.csv')
data = data[(data['strain']=='dilution') & (data['repressors'] > 0) & 
            (data['fold_change'] >= 0)]

# Group the data by each date, run number, replicate, and ATC to compute the means. 
grouped = data.groupby(['carbon', 'temp', 'date', 
                        'run_number', 'atc_ngml']).mean().reset_index()

#%%
# Load the stan model 
model = mwc.bayes.StanModel('../stan/DNA_binding_energy.stan', force_compile=True) 

#%%
# Assign the data dictionary and sample the model. 
rename_params = ['epRA', 'sigma']
sample_dfs = []
summary = []
for g, d in tqdm.tqdm(grouped.groupby(['carbon', 'temp'])):
    data_dict = {'N':len(d),
                 'Nns':4.6E6,
                 'repressors':d['repressors'],
                 'foldchange':d['fold_change']}
    fit, samples = model.sample(data_dict, iter=5000, 
                    control=dict(adapt_delta=0.99))
    params = model.summarize_parameters()
    params['carbon'] = g[0]
    params['temp'] = g[1]
    summary.append(params)
    for r in rename_params:
        df = pd.DataFrame([])
        df['value'] = samples[f'{r}']
        df['parameter'] = r
        df['carbon'] = g[0]
        df['temp'] = g[1]
        sample_dfs.append(df)

samples = pd.concat(sample_dfs)
params = pd.concat(summary)

# Save the dataframes to disk. 
params.to_csv('../../data/DNA_binding_energy_summary.csv', index=False)
samples.to_csv('../../data/DNA_binding_energy_samples.csv', index=False)

#%%
params


# %%
