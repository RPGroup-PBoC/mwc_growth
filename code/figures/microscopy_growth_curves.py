# -*- coding: utf-8 -*- 
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.stats
import mwc.io
import mwc.viz
colors = mwc.viz.pub_style()

growth_data = pd.read_csv('../../data/compiled_growth_microscopy.csv')
growth_data = growth_data[growth_data['area'] >0]
growth_samples = pd.read_csv('../../data/mean_area_growth_samples.csv')

color_key = {c:colors[i] for i, c in enumerate(growth_data['carbon'].unique())}

grouped = growth_data.groupby(['carbon', 'time_min']).mean().reset_index()
fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim([0, 480])
ax.set_ylim([0, 10])
ax.set_xlabel('time [min]', fontsize=12)
ax.set_ylabel(r'segmented area [µm$^2$]', fontsize=12)
# ax.set_yscale('log')
# Plot the data
for g, d in grouped.groupby(['carbon']):
    if g == 'acetate': 
        time = np.arange(0, len(d), 1) * 10
    else:
        time = d['time_min'] - d['time_min'].min()
    ax.plot(time, d['area'] / d.iloc[0]['area'], '.', color=color_key[g], lw=2, label=g)
    
    time_range = np.linspace(0, 500, 500)
    r_median = np.exp(np.mean(growth_samples[f'log_r.{g}']))
    fit = np.exp(time_range * r_median)
    ax.plot(time_range, fit, '-', lw=1, color=color_key[g], label='__nolegend__')
                   
        
        

    
    
