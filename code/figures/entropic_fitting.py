#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import mwc.stats
colors, color_list = mwc.viz.bokeh_theme()
mwc.viz.personal_style()

# Load the foldchange and inference data
foldchange = pd.read_csv('../../data/analyzed_foldchange.csv')
samples = pd.read_csv('../../data/entropic_parameter_samples.csv')

# Restrict the fold-change data
foldchange = foldchange[(foldchange['strain'] == 'dilution') & 
                        (foldchange['carbon'] == 'glucose') & 
                        (foldchange['repressors'] > 10)].copy()

# Define the range of repressors. 
rep_range = np.logspace(1, 4, 100)

# Define the reference parameters
ref_epRA = -13.9 # in kT
ref_epAI = 4.5 # in kT
ref_temp = 37 + 273.15 # in K

#%%
# Set up  the figure canvas
fig, ax = plt.subplots(3, 3, figsize=(5,5), dpi=100)

for a in ax.ravel():
    a.set_yscale('log')
    a.set_xscale('log')
    a.set_xlim([10, 3E3])
    a.set_ylim([1E-3, 1.1])

for i in range(2):
    for j in range(3):
        ax[i, j].xaxis.set_ticklabels([])
        ax[j, i + 1].yaxis.set_ticklabels([])

for i in range(3):
    ax[i, 0].set_ylabel('fold-change', fontsize=8)
    ax[-1, i].set_xlabel('repressors / cell', fontsize=8)

# Assign the temperatures
temp_key = {32:0, 37:1, 42: 2}
edge_colors = {32:colors['dark_blue'], 37:colors['dark_purple'], 
              42:colors['dark_red']}
fill_colors = {32:colors['light_blue'], 37:colors['light_purple'],
               42: colors['light_red']}

for k, v in temp_key.items():
    ax[0, v].set_title(f'{k}°C - prediction', fontsize=10, loc='left', style='italic')

    ax[v, 0].text(-0.5, 0.1, f'{k}°C - inference', rotation='vertical', fontsize=10, style='italic', 
            transform=ax[v, 0].transAxes)

for k_fit, v_fit in temp_key.items():
    for k_pred, v_pred in temp_key.items():
        _ax = ax[v_fit, v_pred]
        # Restrict the fold-change data to the correct temperature. 
        d = foldchange[foldchange['temp']==k_pred]

        # Plot the single cell measurements
        _ax.plot(d['repressors'], d['fold_change'], '.', ms=0.25, 
                color='#c2c2c2', alpha=0.75, rasterized=True)

        # Group by the ATC concentration and compute the mean and SEM
        _g = d.groupby(['date', 'run_number', 'atc_ngml']).mean().reset_index()
        grouped = _g.groupby(['atc_ngml']).agg(('mean', 'sem')).reset_index()

        # Determine the colors
        if k_fit == k_pred:
            fc = 'white'
            _ax.set_facecolor('#ffffff')
        else:
            fc = fill_colors[k_pred]
        ec = edge_colors[k_pred] 

        # Plot the summarized data
        _ax.errorbar(grouped['repressors']['mean'], grouped['fold_change']['mean'],
                      xerr=grouped['repressors']['sem'], 
                      yerr=grouped['fold_change']['sem'],
                      fmt='.', capsize=2, lw=0.75, linestyle='none', color=ec, 
                      markerfacecolor=fc, markeredgecolor=ec, markeredgewidth=0.5,
                      zorder=100)

        # Compute the rescaled temperature adn relative temp
        temp_K = k_pred + 273.15
        _samples = samples[samples['temp']==k_fit]

        # Compute the corrected energies. 
        epRA_star =  ref_epRA -\
             temp_K * _samples[_samples['parameter']=='delta_S_DNA']['value'].values
        epAI_star = rel_temp * ref_epAI -\
             temp_K * _samples[_samples['parameter']=='delta_S_ALLO']['value'].values

        # Define a vector to plot the credible region
        cred_region = np.zeros((2, len(rep_range)))
        for i, r in enumerate(rep_range):
            pact = (1 + np.exp(-epAI_star))**-1
            fc = (1 + pact * (r / 4.6E6) * np.exp(-epRA_star))**-1
            cred_region[:, i] = mwc.stats.compute_hpd(fc, mass_frac=0.95)

        # Plot a shaded band for the credible region
        _ax.fill_between(rep_range, cred_region[0, :], cred_region[1, :],
                        color=fill_colors[k_fit], alpha=0.75, zorder=99)
plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig('../../figs/entropic_parameter_pairwise_prediction.pdf',
            bbox_inches='tight', facecolor='w')


#%%
