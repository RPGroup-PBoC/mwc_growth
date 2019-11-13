#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
colors, _ = mwc.viz.personal_style()


#%% Load the fluctuations and the snap infos
flucts = pd.read_csv('../../data/analyzed_fluctuations.csv')
snaps = pd.read_csv('../../data/analyzed_foldchange.csv')
snaps = snaps[snaps['atc_ngml']==10]
condition_colors = {'glucose':'purple', 'glycerol':'green', 'acetate':'brown',
                    42:'red', 32:'blue'}

#%%
bins = np.linspace(0, 10, 50)
fig, ax = plt.subplots(2, 5, figsize=(7, 3), dpi=100)

iter = 0
for g, d in snaps.groupby(['carbon', 'temp']):
    # Get the corresponding fluctuations
    _flucts = flucts[(flucts['carbon']==g[0]) & (flucts['temp']==g[1])]

    if g[1] == 37:
        tc = colors[condition_colors[g[0]]]
        fc = colors[f'pale_{condition_colors[g[0]]}']
    else:
        tc = colors[condition_colors[g[1]]]
        fc = colors[f'pale_{condition_colors[g[1]]}']
    mwc.viz.titlebox(ax[0,iter], text=f'{g[0].upper()}, {g[1]}° C',
                     size=6, color=tc, bgcolor=fc, boxsize="15%")

    # plot the histogram of the data lengths. 
    _ = ax[0, iter].hist(d['length_um'], bins=bins, histtype='stepfilled',
                        edgecolor=tc, facecolor=fc, alpha=0.75, density=True)
    _ = ax[0, iter].hist(_flucts[['length_1_birth', 'length_2_birth']].values.flatten(), 
                        bins=bins, histtype='stepfilled',
                        edgecolor=colors['black'], facecolor=colors['black'], 
                        alpha=0.25, density=True)
    # Plot the ecdf. 

    # Plot the ecdf. 
    x, y = np.zeros(len(d) + 2), np.zeros(len(d) + 2)
    x[-1] = 100
    y[-1] = 1
    _x, _y = np.sort(d['length_um']), np.arange(1, len(d)+1, 1) / len(d)
    x[1:-1] = _x
    y[1:-1] = _y

    _ = ax[1, iter].step(x, y, color=tc)
    iter += 1

for a in ax.ravel():
    a.set_xlim([0, 6])
for i in range(5):
    ax[0, i].set_xticklabels([])
    ax[0, i].set_ylim([0, 0.7])
    if i > 0:
        for j in range(2):
            ax[j, i].set_yticklabels([])

plt.subplots_adjust(wspace=0.1, hspace=0.1)

# %%
