# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import mwc.stats
import mwc.model
import mwc.viz
colors = mwc.viz.pub_style()

# Load the fold-change data
fc_data = pd.read_csv('../../data/compiled_fold_change.csv')

# Assign colors. 
color_key = {'acetate':colors[1], 'glycerol':colors[2], 'glucose': colors[0]}

# Compute the mean and sem for each carbon. 
dilution = fc_data[(fc_data['strain']=='dilution') & (fc_data['repressors'] >= 20) & (fc_data['fold_change'])]
grouped = dilution.groupby(['carbon', 'atc_ngml', 'date', 'run_number']).mean().reset_index()

# Compute the theoretical prediction. 
rep_range = np.logspace(0, 4, 200)
theo = mwc.model.SimpleRepression(rep_range, ep_r=-13.9, ep_ai=4.5, ka=139, ki=0.53,
                                 effector_conc=0).fold_change()

# Set up the figure canvas and format. 
fig, ax = plt.subplots(1, 1, figsize=(9, 6))
ax.set_xlabel('repressors per cell', fontsize=18)
ax.set_ylabel('fold-change', fontsize=18)
ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([10, 1E3])

# Plot the repressors. 
_ = ax.plot(rep_range, theo, 'k-', lw=2, label='prediction')

# Plot the data.
for g, d in grouped.groupby(['carbon']):
    # Compute the mean and sem for each concentration 
    mean_vals = d.groupby(['atc_ngml'])[['fold_change', 'repressors']].mean()
    sem_vals = d.groupby(['atc_ngml'])[['fold_change', 'repressors']].sem()
                            
    # Plot the points and errors. 
    _ = ax.errorbar(mean_vals['repressors'], mean_vals['fold_change'], 
                   yerr=sem_vals['fold_change'], xerr=sem_vals['repressors'],
                   fmt='o', lw=1, capsize=2, color=color_key[g], ms=5, label=g) 
_ = ax.legend(fontsize=14, loc='lower left')
plt.savefig('../../figs/carbon_repressor_titration.svg', bbox_inches='tight')