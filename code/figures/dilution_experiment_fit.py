#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import phd.viz
import mwc.stats
colors, color_list = mwc.viz.bokeh_theme()


# Load the lineages
data = pd.read_csv('../../data/analyzed_lineages_proper.csv')
data['summed'] = data['I_1_sub'] + data['I_2_sub']
data['fluct'] = (data['I_1_sub'] - data['I_2_sub'])**2


#%%
data.head()



#%%

# Instantiate the figure
fig, ax = plt.subplots(1, 1, figsize=(4, 3))
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('$I_1 + I_2$ [a.u.]', fontsize=10, style='italic')
ax.set_ylabel(r'$\left(I_1 - I_2\right)^2$ [a.u.]', fontsize=10, style='italic')

# Plot the glucose divisions and rasterize
gluc_37 = data[(data['carbon']=='acetate')]


ax.plot(gluc_37['summed'], gluc_37['fluct'], '.', ms=1, alpha=0.25, color='k',
        label='single division', rasterized=True)

bins = 75
binned = mwc.stats.bin_by_events(gluc_37, bins)
ax.errorbar(binned['summed'], binned['fluct'], xerr=binned['summed_sem'], 
           yerr=binned['fluct_sem'], fmt='.', linestyle='none',
           markerfacecolor=colors['light_purple'], markeredgecolor=colors['purple'],
           ms=8, label='75 divisions / bin', color=colors['purple'], capsize=1.5,
           alpha=0.75)

I_tot_range = np.logspace(1, 5.8)
ax.fill_between(I_tot_range, gluc_37['alpha_min'] * I_tot_range, 
            gluc_37['alpha_max'] * I_tot_range, color=colors['purple'], 
       label='$\alpha = %.0f$ [a.u.] / molecule' %gluc_37['alpha_mean'].median())
#%%


#%%
