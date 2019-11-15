#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import mwc.viz
import mwc.model
import mwc.stats
colors, _ = mwc.viz.personal_style()

#%


#%
# Load the data sets and restrict to the carbon sources
data = pd.read_csv('../../data/analyzed_foldchange.csv')
samps = pd.read_csv('../../data/entropic_parameter_samples.csv')

data = data[(data['strain']=='dilution') & (data['repressors'] > 0) & 
            (data['fold_change'] >= 0) & (data['carbon'] == 'glucose') & 
            (data['temp'] != 37)]

summary = data.groupby(['temp', 'date', 'run_number', 'atc_ngml']).mean().reset_index()
summary = summary.groupby(['temp', 'atc_ngml']).agg(('mean', 'sem')).reset_index()

samps = samps[(samps['carbon']=='glucose') & (samps['temp'] != 37)]

# Define the constants for plotting
rep_range = np.logspace(0, 3, 100)

# %%
# Set up the figure canvas
fig, ax = plt.subplots(2, 2, figsize=(3, 3), sharex=True, sharey=True, 
                        dpi=100)
for a in ax.ravel():
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlim([1, 800])
    a.set_ylim([1E-3, 1])

for i in range(2):
    ax[-1, i].set_xlabel('repressors per cell')
    ax[i, 0].set_ylabel('fold-change')

for i in range(2):
    ax[i, i].set_facecolor('white')
    ax[i, i].grid(color=colors['grey'], lw=1)

for a in ax.ravel():
    a.set_xlim([0.9, 1000])
titles = ['32 °C', '42 °C']
temps = [32, 42]
ref_temp = 37 + 273.15
title_colors = [colors['dark_blue'], colors['dark_red']]
bg_colors = [colors['pale_blue'], colors['pale_red']]
face_colors = [colors['light_blue'], colors['light_red']]
for i in range(2):
    mwc.viz.ylabelbox(ax[i, 0], titles[i], title_colors[i], bgcolor=bg_colors[i],
                        size=6, boxsize="15%")
    mwc.viz.titlebox(ax[0, i], titles[i], title_colors[i], bgcolor=bg_colors[i],
                        size=6, boxsize="15%")
    
    if i > 0:
        # apply offset transform to all y ticklabels.
        dx = -13 / fig.dpi
        dy = 0
        offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
        for label in ax[i, 0].yaxis.get_majorticklabels():
            label.set_transform(label.get_transform() + offset)


# Plot the predictions
for i, fit_temp in enumerate(temps):
    # Compute the credible region for each prediction 
    true_epRA = samps[(samps['parameter']=='true_epRA') & (samps['temp']==fit_temp)]['value'].values
    true_epAI = samps[(samps['parameter']=='true_epAI') & (samps['temp']==fit_temp)]['value'].values
    del_S = samps[(samps['parameter']=='delta_S') & (samps['temp']==fit_temp)]['value'].values
    allo_S = samps[(samps['parameter']=='delta_S_vib') & (samps['temp']==fit_temp)]['value'].values

    for j, pred_temp in enumerate(temps):
        # if i == j:
            # epRA = samps[(samps['parameter']=='epRA_star') & (samps['temp']==fit_temp)]['value'].values
            # epAI = samps[(samps['parameter']=='epAI_star') & (samps['temp']==fit_temp)]['value'].values
        # else:
        pred_temp += 273.15
        epRA = (ref_temp/ pred_temp) * true_epRA  - pred_temp * del_S
        epAI =  (ref_temp / pred_temp ) * true_epAI  - pred_temp * allo_S

        cred_region = np.zeros((2, len(rep_range)))
        pref = np.exp(-epRA) / (1 + np.exp(-epAI))
        for k, r in enumerate(rep_range):
            theo = (1 + (r / 4.6E6) * pref)**-1
            cred_region[:, k] = mwc.stats.compute_hpd(theo, 0.95)

        ax[i, j].fill_between(rep_range, cred_region[0, :], cred_region[1, :], 
                color=title_colors[i], alpha=0.25)

# Plot the data
for i, temp1 in enumerate(temps):
    for j, temp2 in enumerate(temps):
        if i == j:
            fill  = 'white'
        else:
            fill = face_colors[i]
        # Isolate the data. 
        d = summary[summary['temp']==temps[i]]
        ax[j, i].errorbar(d['repressors']['mean'], d['fold_change']['mean'],
                        xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'],
                        color=title_colors[i], fmt='.', ms=8, markerfacecolor=fill,
                        markeredgewidth=0.75, linestyle='none', capsize=1,
                        lw=0.75)
plt.subplots_adjust(wspace=0.05, hspace=0.05)

# %%


# %%
