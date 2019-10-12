#%%
import numpy as np
import pandas as pd
import mwc.bayes
import mwc.stats
import tqdm

# Load the data
data = pd.read_csv('../../data/analyzed_foldchange.csv')

# Isolate the data to the "true" strains and compute the summary statistics
data = data[(data['strain']=='dilution') & (data['repressors'] >= 10) & (data['fold_change'] > -0.001)]

# Compute the replicate summary statistics
replicate = data.groupby(['date', 'run_number', 
                          'carbon', 'temp', 'atc_ngml']).mean().reset_index()

# Load the statistical model and compile if need be
model = mwc.bayes.StanModel('../stan/empirical_F_inference.stan') 
                            # force_compile=True)
#%%
# Instantiate storage vectors and begin the inference for each condition
summary_dfs = []
for g, d in  tqdm.tqdm(replicate.groupby(['carbon', 'temp'])):
    d = d.copy()

    # Generate the identification vector
    d['idx'] = d.groupby('atc_ngml').ngroup() + 1   

    # Create the mapping from idx number to atc concentration.
    _d = d.groupby('idx')['atc_ngml'].mean().reset_index()
    atc_map = {i:a for i, a in zip(_d['idx'].values, _d['atc_ngml'].values)}

    # Set up the data dictionary. 
    data_dict = {'N':len(d), 'J':d['idx'].max(), 
                 'idx':d['idx'], 'foldchange':d['fold_change'],
                 'repressors':d['repressors']}
    fit, samples = model.sample(data_dict, iter=5000, 
                                control=(dict(adapt_delta=0.99)))

    # Summarize the parameters
    params = model.summarize_parameters()
    for k, v in atc_map.items():
        params.loc[params['dimension']==k, 'atc_ngml'] = v

    params['carbon'] = g[0]
    params['temp'] = g[1]

    # add the parameter info to the summary list
    summary_dfs.append(params)

# Concatenate the summaries into a single tidy dataframe.
df = pd.concat(summary_dfs)
df.to_csv('../../data/inferred_empirical_F.csv', index=False)
#%%
