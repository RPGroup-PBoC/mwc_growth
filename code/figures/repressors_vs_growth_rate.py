# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import mwc.stats
import bokeh.palettes
mwc.viz.pub_style()
colors = bokeh.palettes.viridis(4)


# Load the various data
fc_data = pd.read_csv('../../data/compiled_fold_change.csv')
growth_samples = pd.read_csv('../../data/mean_area_growth_samples.csv')

# Determine the HPD of growth rate.  ind = np.argmax(growth_samples['lp__'].values)
growth_rates = {}
for i, c in enumerate(fc_data['carbon'].unique()):
    r_mode = growth_samples.iloc[ind][f'r.{c}']
    r_min, r_max = mwc.stats.compute_hpd(growth_samples[f'r.{c}'], 0.95)
    fc_data.loc[fc_data['carbon']==c, 'rate_mode'] = 1 / r_mode 
    fc_data.loc[fc_data['carbon']==c, 'rate_min'] =  1 / r_min
    fc_data.loc[fc_data['carbon']==c, 'rate_max'] =  1 / r_max
        
# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(4, 5))
ax.set_xlabel('growth rate [min$^{-1}$]', fontsize=18)
ax.set_ylabel('repressors per cell', fontsize=18)
ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)


concs = fc_data[(fc_data['strain']=='dilution') & ((fc_data['atc_ngml']==1) |
                                                   (fc_data['atc_ngml']==3) |
                                                   (fc_data['atc_ngml']==10))]
color_key = {c:colors[i] for i, c in enumerate(concs['atc_ngml'].unique())}
for g, d in concs.groupby(['atc_ngml']): 
    grouped = d.groupby(['carbon', 'date', 'run_number']).mean().reset_index()
    grouped_mean = grouped.groupby(['carbon']).mean().reset_index()
    grouped_sem = grouped.groupby(['carbon']).sem().reset_index()
    grouped_mean = grouped_mean.sort_values('rate_mode')
    _ = ax.errorbar(grouped_mean['rate_mode'], grouped_mean['repressors'], grouped_sem['repressors'],
                    fmt='o', lw=1, color=color_key[g], label=int(g), ms=5, linestyle='-', zorder=100,
                   capsize=2)

leg = ax.legend(title='ATC [ng / mL]', fontsize=12)    
leg.get_title().set_fontsize(12)
plt.savefig('../../figs/reps_v_growthrate.svg', bbox_inches='tight') 