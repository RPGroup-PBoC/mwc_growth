#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.stats
import mwc.viz
colors, color_list = mwc.viz.bokeh_theme()
mwc.viz.personal_style()


# Load the data
growth_stats = pd.read_csv('../../data/compiled_growth_statistics.csv')
epRA_summary = pd.read_csv('../../data/DNA_binding_energy_summary.csv')
entropic_summary = pd.read_csv('../../data/entropic_parameter_summary.csv')

# Restrict the parameters in the epRA summary to only the DNA binding energy. 
epRA_summary = epRA_summary[epRA_summary['parameter']=='epRA'].copy()

# Compute the growth rates and include in the binding energy dataframe
for g, d in growth_stats.groupby(['carbon', 'temp']):
    growth_rate = d[d['parameter']=='max df']['value'].mean()* 60 # in hr^-1
    dbl_time = d[d['parameter']=='inverse max df']['value'].mean() # in min^-1

    # Add the values to the epRA summary dataframe.
    epRA_summary.loc[(epRA_summary['carbon']==g[0]) &
                     (epRA_summary['temp']==g[1]), 'mean_k'] = growth_rate
    epRA_summary.loc[(epRA_summary['carbon']==g[0]) &
                     (epRA_summary['temp']==g[1]), 'mean_dbl'] = dbl_time

#%%
# Set up the figure canvas. 
fig, ax = plt.subplots(1, 1, figsize=(3, 2), dpi=100)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
ax.set_ylabel('inferred DNA binding energy [$k_BT$]', fontsize=8)
ax.set_ylim([-16, -12])

# Define the colors and glyphs
glyphs = {32: 's', 37: 'o', 42: 'v'}
edge_colors = {'glucose': colors['dark_purple'], 
               'glycerol':colors['dark_green'],
               'acetate': colors['dark_orange']}
fill_colors = {'glucose': colors['light_purple'], 
               'glycerol':colors['light_green'],
               'acetate': colors['light_orange']}

# Loop through each condition and plot
for g, d in epRA_summary.groupby(['carbon', 'temp']):
    ax.vlines(d['mean_k'], d['hpd_min'], d['hpd_max'], lw=0.35, 
                color=edge_colors[g[0]], label='__nolegend__')
    ax.plot(d['mean_k'], d['median'], marker=glyphs[g[1]], 
            markeredgecolor= edge_colors[g[0]], markerfacecolor=fill_colors[g[0]],
            markeredgewidth=0.25,  label=f'{g[0]}, {int(g[1])}° C', linestyle='none',
            ms=4)


# Fill between the "true" fit value
ax.fill_between(np.linspace(0.1, 0.7), -14.1, -13.7, color=colors['light_grey'], alpha=0.5,
zorder=1, label='Garcia & Phillips, 2011')
ax.legend(ncol=3, bbox_to_anchor=(1.1, 1.3))
plt.savefig('../../figs/binding_energy_V_growthrate.pdf', bbox_inches='tight',
           facecolor='white')


# %%
# Make a figure of inferred binding energy vs temperature. 
fig, ax = plt.subplots(1, 1, figsize=(3, 2), dpi=100)

#Look only at the temperature. 
temp_data = epRA_summary[epRA_summary['carbon']=='glucose']
temp_range = np.linspace(300, 320, 100)
epRA_ref = -13.9
ls = {32: ':', 37: '-', 42: '--'}
for g, d in temp_data.groupby(['temp', 'median']):
    ax.plot(1000 / (g[0] + 273.15), g[1], marker=glyphs[g[0]], 
            markerfacecolor=colors['light_purple'],
            markeredgecolor=colors['dark_purple'],
            linestyle='None', label=f'{g[0]}° C',
            markeredgewidth=0.25, ms=4)
    ax.vlines(1000 / (g[0] + 273.15), d['hpd_min'], d['hpd_max'], 
              color=colors['dark_purple'], lw=0.35, label='__nolegend__')

    # Get whatever teh entropic component is. 
    dS = entropic_summary[(entropic_summary['temp']==g[0])  & 
                        (entropic_summary['parameter']=='delta_S_DNA')]
    # ax.plot(temp_range, epRA_ref - temp_range * dS['median'].values[0], 
    #         linestyle=ls[g[0]], color=colors['dark_purple'],
    #         label=r'$\Delta\varepsilon_{RA}^{(37^\circ C)} - T\Delta S^{(%s^\circ C)}$' % int(g[0])) 

ax.legend(ncol=3, bbox_to_anchor=(1.1, 1.4))
ax.set_xlabel('temperature [K]')
ax.set_ylabel('inferred DNA binding energy [$k_BT$]')
plt.savefig('../../figs/inferred_binding_energy_v_T.pdf', 
bbox_inches='tight', facecolor='white')
#%%
