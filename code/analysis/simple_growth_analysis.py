# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import pandas as pd
import numpy as np
import mwc.io
import mwc.stats
import mwc.bayes

# Load the growth curve data. 
growth = pd.read_csv('../../data/compiled_growth_microscopy.csv')
grouped = growth.groupby(['carbon', 'time_min']).mean().reset_index()

# Add proper identifiers. 
grouped['idx'] = grouped.groupby(['carbon']).ngroup() + 1
dfs = []
for g, d in grouped.groupby(['carbon']):
    d = d.copy()
    d.sort_values(by='time_min', inplace=True)
    d['time_min'] -= d['time_min'].min()
    if g == 'acetate':
        d['time_min'] = np.arange(0, len(d), 1) * 10
    dfs.append(d)
grouped = pd.concat(dfs)

grouped['area'] *= 0.065**2

# Assemble the data dictionary. 
data_dict = {'J':np.max(grouped['idx']),
            'N': len(grouped),
            'area': grouped['area'],
            'time': grouped['time_min'],
            'index_1': grouped['idx']}

model = mwc.bayes.StanModel('../stan/projected_mean_growth.stan', 
                            data_dict=data_dict, force_compile=True)
samples = model.sample()
samples, samples_df = samples

label_key = {grouped[grouped['carbon']==c]['idx'].unique()[0]:c for c in grouped['carbon'].unique()}
new_names = {k:f"{k.split('[')[0]}.{label_key[int(float(k.split('[')[1][0]))]}" for k in samples.flatnames}
samples_df.rename(columns=new_names, inplace=True)
samples_df.to_csv('../../data/mean_area_growth_samples.csv', index=False)

