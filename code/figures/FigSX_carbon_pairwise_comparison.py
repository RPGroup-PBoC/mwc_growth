#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import mwc.viz
import mwc.model
import mwc.stats
colors, _ = mwc.viz.personal_style()

# %%
# Load the data sets and restrict to the carbon sources
data = pd.read_csv('../../data/analyzed_foldchange.csv')
stats = pd.read_csv('../../data/DNA_binding_energy_summary.csv')

data = data[(data['strain']=='dilution') & (data['repressors'] > 0) & 
            (data['fold_change'] >= 0) & (data['temp'] == 37)]
summary = data.groupby(['carbon', 'date', 'run_number', 'atc_ngml']).mean().reset_index()
summary = summary.groupby(['carbon', 'atc_ngml']).agg(('mean', 'sem')).reset_index()
stats = stats[(stats['temp']==37)]

# Define the constants for plotting
rep_range = np.logspace(0, 3, 100)
# %%
# Set up the figure canvas
fig, ax = plt.subplots(3, 3, figsize=(5.5, 5.5), sharex=True, sharey=True, 
                        dpi=100)
for a in ax.ravel():
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlim([1, 800])

for i in range(3):
    ax[-1, i].set_xlabel('repressors per cell')
    ax[i, 0].set_ylabel('fold-change')

for i in range(3):
    ax[i, i].set_facecolor('white')
    ax[i, i].grid(color=colors['grey'], lw=1)

titles = ['acetate', 'glycerol', 'glucose']
title_colors = [colors['dark_brown'], colors['dark_green'], colors['dark_purple']]
bgcolors = [colors['pale_brown'], colors['pale_green'], colors['pale_purple']]
face_colors = [colors['brown'], colors['light_green'], colors['light_purple']]
for i in range(3):
    mwc.viz.ylabelbox(ax[i, 0], titles[i].upper(), title_colors[i], bgcolor=bgcolors[i],
                size=6, boxsize="15%")
    mwc.viz.titlebox(ax[0, i], titles[i].upper(), title_colors[i], bgcolor=bgcolors[i],
                size=6, boxsize="15%")
    
    if i > 0:
        # apply offset transform to all y ticklabels.
        dx = -13 / fig.dpi
        dy = 0
        offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
        for label in ax[i, 0].yaxis.get_majorticklabels():
            label.set_transform(label.get_transform() + offset)


# Plot the predictions
for i, pred in enumerate(titles):
    # Get the binding energy values for the prediction strain 
    low, high = stats[(stats['carbon']==pred) & 
                      (stats['parameter']=='epRA')][
                      ['hpd_min', 'hpd_max']].values[0]
    # Compute the theory
    theo_min = mwc.model.SimpleRepression(R=rep_range, ep_r=low, ka=139, ki=0.53,
                                          ep_ai=1000, effector_conc=0).fold_change()
    theo_max = mwc.model.SimpleRepression(R=rep_range, ep_r=high, ka=139, ki=0.53,
                                          ep_ai=1000, effector_conc=0).fold_change()
    for j, fit in enumerate(titles): 
        ax[i, j].fill_between(rep_range, theo_min, theo_max, color=title_colors[i],
                alpha=0.25)

# Plot the data
for i, carb in enumerate(titles):
    for j in range(3):
        if i == j:
            fill  = 'white'
        else:
            fill = face_colors[i]
        # Isolate the data. 
        d = summary[summary['carbon']==carb]
        ax[j, i].errorbar(d['repressors']['mean'], d['fold_change']['mean'],
                        xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'],
                        color=title_colors[i], fmt='.', ms=8, markerfacecolor=fill,
                        markeredgewidth=0.75, linestyle='none', capsize=1,
                        lw=0.75)
plt.subplots_adjust(wspace=0.05, hspace=0.05)
plt.savefig('../../figs/FigS_carbon_binding_energy_pairwise_fc.pdf', 
            bbox_inches='tight', facecolor='white')

# %%
