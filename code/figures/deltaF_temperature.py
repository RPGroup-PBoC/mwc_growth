#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import mwc.stats
colors, color_list = mwc.viz.personal_style()

# Define the repressor copy number range. 
rep_range = np.logspace(0, 4, 200)
ref_rep = 100

# Define thermodynamic constants
ep_ai = 4.5 # in kT
ep_RA = -13.9 # in kT
T_ref = 37 + 273.15 # in K
Nns = 4.6E6 # in bp
sigma = 0.2 # Error in estimation of binding energy

# Load the actual data. 
fc_data = pd.read_csv('../../data/analyzed_foldchange.csv')
inferred_F = pd.read_csv('../../data/inferred_empirical_F.csv')
fc_data = fc_data[(fc_data['fold_change'] > -0.001) & (fc_data['carbon']=='glucose')]

# Isolate the fc data to the relevant measurements 
fc_data = fc_data[(fc_data['strain']=='dilution') & 
                  (fc_data['repressors'] >= 10)].copy()
inferred_F = inferred_F[inferred_F['carbon']=='glucose'].copy()

# Compute the summary statistics
rep_summary = fc_data.groupby(['date', 'run_number', 
                               'atc_ngml', 'temp']).mean().reset_index()
summary = rep_summary.groupby(['atc_ngml', 'temp']).agg(('mean', 'sem')).reset_index()


# Compute the reference bits
ref_pact = (1 + np.exp(-ep_ai))**-1
ref_fc = (1 + ref_pact * rep_range * np.exp(-ep_RA) / Nns)**-1

# Define the temperature colors. 
colors_fill = {32:colors['light_blue'], 37:colors['light_purple'], 42:colors['light_red']}
colors_edge = {32:colors['dark_blue'], 37:colors['dark_purple'], 42:colors['dark_red']}

#%%
# Set up the figure canvas
fig, ax = plt.subplots(2, 3, figsize=(4.5, 3), dpi=150)
for i in range(3):
    ax[0, i].set_xscale('log')
    ax[0, i].set_yscale('log')
    ax[0, i].set_xlim([10, 700])
    ax[0, i].set_ylim([1E-3, 1])


# Define the axes for temperatures.
temp_axes = {32:0, 37:1, 42:2}

# Plot the theory. 
for T, i in temp_axes.items():
    rel_T = (T + 273.15) / T_ref

    # Compute the theoretical fold-change
    pact = (1 + np.exp(-rel_T * ep_ai))**-1
    fc = (1 + pact * rep_range * np.exp(-rel_T * ep_RA) / Nns)**-1

    # Compute the deltaF 
    delF = -np.log(pact / ref_pact) - np.log(rep_range / ref_rep) + (ep_RA * (1 - rel_T))
    if i != 1:
        ms =  '--'
    else:
        ms = '-'
    ax[0, i].plot(rep_range, fc, 'k', linestyle=ms, linewidth=0.75)
    ax[1, i].plot(rep_range / ref_rep, -delF, 'k', linestyle=ms, linewidth=0.75)

ax[1, 0].set_xlim([0.1, 2])
ax[1, 1].set_xlim([0.1, 3])
ax[1, 2].set_xlim([0.1, 7])


# Plot the data
for g, d in summary.groupby(['temp']):
    ax[0, temp_axes[g]].errorbar(d['repressors']['mean'], d['fold_change']['mean'], 
            xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'], fmt='.', 
            linestyle='none', markersize=8, markerfacecolor=colors_fill[g],
            markeredgecolor=colors_edge[g], markeredgewidth=0.75, 
            color=colors_fill[g])


# Plot the empirical delta F vs the repressor copy number
for g, d in inferred_F.groupby(['temp']):
    _reps = summary[summary['temp']==g]['repressors']
    F_median = d[d['parameter']=='empirical_F']['median'].values
    F_min = d[d['parameter']=='empirical_F']['hpd_min'].values
    F_max = d[d['parameter']=='empirical_F']['hpd_max'].values
    F_theo = -np.log(ref_pact) - np.log(ref_rep / Nns) + ep_RA
    delF_median = F_median - F_theo
    delF_max = F_max - F_theo
    delF_min = F_min - F_theo

    ax[1, temp_axes[g]].vlines(_reps['mean'] / ref_rep, -delF_min, -delF_max, linewidth=0.75,
                            color=colors_edge[g])
    ax[1, temp_axes[g]].errorbar(_reps['mean'] / ref_rep, -delF_median, xerr=_reps['sem'] / ref_rep,
                linestyle='none', color=colors_edge[g], markerfacecolor=colors_fill[g],
                markeredgecolor=colors_edge[g], markeredgewidth=0.75, ms=8, fmt='.')

#%%
