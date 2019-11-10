#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import mwc.model
import mwc.stats
colors, color_list = mwc.viz.personal_style()

# Load the data sets and restrict
samples = pd.read_csv('../../data/DNA_binding_energy_samples.csv')
stats = pd.read_csv('../../data/DNA_binding_energy_summary.csv')
samples = samples[samples['temp'] == 37]
stats = stats[stats['temp'] == 37]

# %%
# Set up the figure canvas
fig, ax = plt.subplots(2,2, figsize=(4, 4), dpi=100)

# Format the axes
ax[0, 1].axis(False)
ax[0, 0].yaxis.grid(False)
ax[0, 0].set_yticks([])
ax[1, 1].yaxis.grid(False)
ax[1, 1].set_yticks([])
ax[0, 0].set_xticklabels([])

# Set the titles of the marginals for clarity
mwc.viz.titlebox(ax[0, 0], '$\epsilon$', boxsize="15%", pad=0.03, color=colors['black'])
mwc.viz.titlebox(ax[1, 1], '$\sigma$', boxsize="15%", pad=0.03, color=colors['black'])

# Add the axis labels
ax[1, 0].set_ylabel('$\sigma$')
ax[1, 0].set_xlabel('$\epsilon$ [$k_BT$]')

carb_colors = {'glucose':colors['purple'], 'glycerol':colors['green'], 
          'acetate':colors['brown']}
for g, d in samples.groupby('carbon'):
    ep = d[d['parameter']=='epRA']['value'].values
    sig = d[d['parameter']=='sigma']['value'].values
    ax[1, 0].plot(ep, sig, ',', color=carb_colors[g], alpha=0.5)

plt.subplots_adjust(wspace=0.05, hspace=0.05)
# %%

samples['parameter'].unique()


# %%
