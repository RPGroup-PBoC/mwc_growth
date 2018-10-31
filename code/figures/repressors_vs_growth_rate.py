# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import mwc.stats
colors = mwc.viz.pub_style()

# Load the various data
fc_data = pd.read_csv('../../data/compiled_fold_change.csv')
growth_samples = pd.read_csv('../../data/mean_area_growth_samples.csv')

# Determine the HPD of growth rate. 
ind = np.argmax(growth_samples['lp__'].values)
growth_rates = {}
for i, c in enumerate(fc_data['carbon']):
    r_mode = growth_samples.iloc[ind][f'r.{c}']
    r_min, r_max = mwc.stats.compute_hpd(growth_samples[f'r.{c}'], 0.95)
    fc_data.loc[fc_data['carbon']==c, 'rate_mode'] = r_mode 
    fc_data.loc[fc_data['carbon']==c, 'rate_min'] = r_min
    fc_data.loc[fc_data['carbon']==c, 'rate_max'] = r_max
        
# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(6, 4))

for g, d in fc_data[fc_data['strain']=='dilution'].groupby(['atc_ngml']):
    grouped = d.groupby(['carbon', 'date', 'run_number']).mean().reset_index()
    
    
    