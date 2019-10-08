#%%
import numpy as np
import pandas as pd
import mwc.bayes
import mwc.stats

# Load the fold-change data and growth rate stats
foldchange = pd.read_csv('../../data/analyzed_foldchange.csv')
stats = pd.read_csv('../../data/compiled_growth_statistics.csv')

# Determine doubling times
stats = stats[((stats['carbon']=='glucose') | (stats['carbon']=='acetate') |
                (stats['carbon']=='glycerol')) & 
                ((stats['temp']==37) |  (stats['temp']==32) | 
                (stats['temp']==42))] 

tidy_stats = pd.DataFrame([])
for g, d in stats.groupby(['date', 'carbon', 'temp', 'run_number']):
    growth_rate = d[d['parameter']=='max df']['value'].values[0]
    growth_err = d[d['parameter']=='max df std']['value'].values[0]
    dbl_time = d[d['parameter']=='inverse max df']['value'].values[0]
    dbl_err = d[d['parameter']=='inverse max df std']['value'].values[0]

    tidy_stats = tidy_stats.append({'date':g[0], 'carbon':g[1], 'temp_C':g[2], 'run_number':g[3],
                                    'growth_rate':growth_rate, 
                                    'dbl_time':dbl_time,
                                    'growth_err':growth_err,
                                    'dbl_err':dbl_err}, 
                                    ignore_index=True)
tidy_stats['growth_rate'] *= 60
tidy_stats['growth_err'] *= 60

# Summarize the growth rates
tidy_stats = tidy_stats.groupby(['carbon', 'temp_C']).agg(('mean', 'sem')).reset_index()
for g, d in tidy_stats.groupby(['carbon', 'temp_C']):
    foldchange.loc[(foldchange['carbon']==g[0]) & (foldchange['temp']==g[1]),
                   'rate_mean'] = d['growth_rate']['mean'].values[0]
    foldchange.loc[(foldchange['carbon']==g[0]) & (foldchange['temp']==g[1]),
                   'rate_sem'] = d['growth_rate']['sem'].values[0]
    foldchange.loc[(foldchange['carbon']==g[0]) & (foldchange['temp']==g[1]),
                   'dbl_mean'] = d['dbl_time']['mean'].values[0]
    foldchange.loc[(foldchange['carbon']==g[0]) & (foldchange['temp']==g[1]),
                   'dbl_sem'] = d['dbl_time']['sem'].values[0]


# Keep only the dilution strains from the fold-change and repressor counts > 0
data = foldchange[(foldchange['strain']=='dilution')  & 
                   (foldchange['repressors'] > 0)].copy()


# %%
# Load and compile the stan model (if needed)
model = mwc.bayes.StanModel('../stan/repressor_vs_growthrate.stan', force_compile=True)

#%%
# Go through both the temp and carbon source variations
carbon = data[data['temp']==37]
temp = data[data['carbon']=='glucose']

# For each, group by the date and run number and compute the mean
carbon_grp = carbon.groupby(['date', 'run_number', 'temp',
                            'carbon', 'atc_ngml', 'rate_mean'])['repressors'].mean().reset_index()
temp_grp = temp.groupby(['date', 'run_number', 'carbon',
                         'temp', 'atc_ngml', 'rate_mean'])['repressors'].mean().reset_index()

# Assign identifiers for each unique ATC concentration
carbon_grp['idx'] = carbon_grp.groupby('atc_ngml').ngroup() + 1
temp_grp['idx'] = temp_grp.groupby('atc_ngml').ngroup() + 1

# Combine the grouped dataframes for iteration
grp = [carbon_grp, temp_grp]
variable = ['carbon', 'temp']
samples_dfs, summary_dfs  = [], []
for i, g in enumerate(grp):
    # Generate the data dictionary
    data_dict = {'J':g['idx'].max(), 
                 'N':len(g), 
                 'idx':g['idx'],
                 'growth_rate':g['rate_mean'],
                 'repressors': g['repressors']}

    # Perform the inference
    fit, samples = model.sample(data_dict)
    params = model.summarize_parameters()

    # Add identifying information
    samples['variable'] = variable[i]
    params['variable'] = variable[i]
    
    # Map the dimension of the intercept parameters to the atc concentration
    dims = params['dimension'].unique()

    for d in params['dimension'].unique():
        conc =  g[g['idx']==d]['atc_ngml'].unique()[0]
        params.loc[params['dimension']==d, 'atc_ngml'] = conc
        samples.rename(columns={f'intercept[{d}]': f'intercept[{conc}_ngml]'})

    samples_dfs.append(samples)
    summary_dfs.append(params)

# Concatenate the frames and save to disk
samples_df = pd.concat(samples_dfs)
summary_df = pd.concat(summary_dfs)
samples_df.to_csv('../../data/linear_regression_samples.csv', index=False)
summary_df.to_csv('../../data/linear_regression_summary.csv', index=False)

#%%
