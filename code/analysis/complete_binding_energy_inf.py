#%%
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import mwc.bayes
import mwc.stats
import mwc.model
import matplotlib.transforms
colors, _ = mwc.viz.personal_style()

# %%
# Load the various data sets. 
data = pd.read_csv('../../data/analyzed_foldchange.csv')
data = data[(data['carbon']=='glucose') & (data['temp']==37) & 
            (data['strain']=='dilution') & (data['repressors'] > 0) & 
            (data['fold_change'] >= 0)]

# Summarize the data
large_only = data[data['size']=='large']
medium_only = data[data['size']=='medium']
small_only = data[data['size']=='medium']
data = data.groupby(['date', 'run_number', 'atc_ngml']).mean().reset_index()
large_only = large_only.groupby(['date', 'run_number', 'atc_ngml']).mean().reset_index()
small_only = small_only.groupby(['date', 'run_number', 'atc_ngml']).mean().reset_index()
medium_only = medium_only.groupby(['date', 'run_number', 'atc_ngml']).mean().reset_index()


# Load the garcia and brewster O2 data
old_gods = pd.read_csv('../../data/Garcia2011_Brewster2014.csv', comment='#')
old_gods = old_gods[old_gods['operator']=='O2']
old_gods.rename(columns={'repressor':'repressors'}, inplace=True)

garcia = old_gods[old_gods['author']=='garcia']
brewster = old_gods[old_gods['author']=='brewster']

# Load the Razo-Mejia O2 data. 
ind_data = pd.read_csv('../../data/RazoMejia_2018.csv', comment='#')
ind_data.rename(columns={'fold_change_A':'fold_change'}, inplace=True)
ind_data['repressors'] *= 2
ind_data = ind_data[(ind_data['operator']=='O2') & (ind_data['repressors'] > 0) &
                    (ind_data['IPTG_uM']==0)]
# %%
# Load the stan model
model = mwc.bayes.StanModel('../stan/DNA_binding_energy.stan')

# %%
# Perform the inference. 
summ_dfs = []
data_dict = {'no_correction':data, 'correction':data, 'large_only':large_only, 
            'all_divided':data, 'medium_only':medium_only, 'small_only':small_only, 'garcia':garcia,
            'brewster':brewster, 'razo-mejia':ind_data}
for k, v in data_dict.items():
    # Define the data dictionary. 
    data_dict = {'N':len(v), 'foldchange':v['fold_change'],
                'Nns':4.6E6}
    if k == 'no_correction':
        data_dict['repressors'] = v['raw_repressors']

    elif k == 'all_divided':
        data_dict['repressors'] = v['raw_repressors'].values * 2

    else:
        data_dict['repressors'] = v['repressors']
    fit, samples = model.sample(data_dict)
    params = model.summarize_parameters()
    params['source'] = k 
    summ_dfs.append(params)
stats = pd.concat(summ_dfs)

# Save the sampling summary to disk. 
stats.to_csv('../../data/R_correction_DNA_binding_energy_summary.csv',
            index=False)

# %%
