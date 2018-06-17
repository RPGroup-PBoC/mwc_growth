# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import pystan
import matplotlib.pyplot as plt
import glob
sys.path.insert(0,'../../../')
import mwc.viz
import mwc.bayes
colors = mwc.viz.personal_style()

# load the data.
data = pd.read_csv(glob.glob('output/*.csv')[0])

#%% Compile the pystan model 
print('compiling model...')
model = pystan.StanModel('../../stan/hierarchical_growth_curve_single_condition.stan')
print('finished!')

#%% Set up the data dictionary and sample.
data_dict = {'J':data['colony_idx'].astype(int).max(), 'N':len(data), 
            'colony':data['colony_idx'].astype(int), 'fractional_area':data['fractional_area'],
            'time':data['time_min']}
print('beginning sampling....')
samples = model.sampling(data_dict, iter=5000, chains=2)
print('finished!')
# %% 
sample_df = mwc.bayes.chains_to_dataframe(samples)
stats = mwc.stats.compute_statistics(sample_df)

# %%
stats
# %%
fig, ax = plt.subplots()
ax.plot(data['time_min'], data['fractional_area'], 'k.')

# %%
data