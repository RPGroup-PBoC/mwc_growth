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
growth_data['area'] *= 0.065**2
growth_samples = pd.read_csv('../../data/mean_area_growth_samples.csv')

color_key = {c:colors[i] for i, c in enumerate(growth_data['carbon'].unique())}

grouped = growth_data.groupby(['carbon', 'time_min']).mean().reset_index()
fig, ax = plt.subplots(figsize=(5, 5))
ax.set_xlim([0, 480])
ax.set_ylim([0, 60])
ax.set_xlabel('time [min]', fontsize=18)
ax.set_ylabel(r'segmented area [µm$^2$]', fontsize=18)
ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)
# ax.set_yscale('log')
# Plot the data
t_double = {}
for g, d in grouped.groupby(['carbon']):
    if g == 'acetate': 
        time = np.arange(0, len(d), 1) * 10
    else:
        time = d['time_min'] - d['time_min'].min()
    ax.plot(time, d['area'] , '.', color=color_key[g], ms=5, label=g)
    
    time_range = np.linspace(0, 500, 500)
    cred_region = np.zeros((2, len(time_range)))
    mode_double = (np.log(2) * growth_samples.iloc[np.argmax(growth_samples['lp__'].values)][f'r.{g}'])
    min_double, max_double = mwc.stats.compute_hpd(np.log(2) * growth_samples[f'r.{g}'], 0.95)
    t_double[g] = [mode_double, mode_double - min_double, max_double - mode_double]
    for i, t in enumerate(time_range):
        fit = np.exp(growth_samples[f'log_A0.{g}']) * np.exp(t / growth_samples[f'r.{g}'])
        cred_region[:, i] = mwc.stats.compute_hpd(fit, 0.95)
    ax.fill_between(time_range, cred_region[0, :], cred_region[1, :], alpha=0.5, color=color_key[g])
                   
plt.savefig('../../figs/microscopy_growth_curves.svg', bbox_inches='tight') 
