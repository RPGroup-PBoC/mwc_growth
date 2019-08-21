#%%
import numpy as np 
import pandas as pd 
import mwc.bayes 
import mwc.stats

snaps = pd.read_csv('../../data/raw_compiled_snaps.csv')

# Choose a single sample. 
snaps['idx'] = snaps.groupby(['date', 'run_number', 'temp', 'carbon']).ngroup() + 1
sel = snaps[snaps['idx']==15]

# Load the model.
model = mwc.bayes.StanModel('../stan/bg_subtraction.stan')

#%%
atc = []
mch = []
auto_list = []
# Set up the data dict
auto = sel[sel['strain']=='auto']
for g, d in sel[sel['strain']=='dilution'].groupby(['atc_ngml']):
    data_dict = {'N':len(auto), 'J':len(d), 'noise':auto['fluor2_mean_death'], 
                'signal':d['fluor2_mean_death']}

    fit, sample = model.sample(data_dict)
    atc.append(g)
    mch.append(np.mean(sample['mu2'].values))
    auto_list.append(np.mean(sample['mu1'].values))

#%%
