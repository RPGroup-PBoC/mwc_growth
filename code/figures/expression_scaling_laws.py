# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import mwc.viz
import seaborn as sns
import joypy
colors, color_list = mwc.viz.bokeh_theme()
mwc.viz.personal_style()

# Load the fold-change data and growth rate stats
foldchange = pd.read_csv('../../data/analyzed_foldchange.csv')
stats = pd.read_csv('../../data/compiled_growth_statistics.csv')


# Define the colors for the conditions
fill_colors = {'acetate': '#e1bb96', 'glycerol': colors['light_green'],
               'glucose':colors['light_purple'], 37: colors['light_purple'],
               32:colors['light_blue'], 42:colors['light_red']}
edge_colors = {'acetate': '#764f2a', 'glycerol': colors['green'],
               'glucose':colors['purple'], 37: colors['purple'],
               32:colors['blue'], 42:colors['red']}

# Assign atc color
sorted_atc = np.sort(foldchange['atc_ngml'].unique())
_colors = sns.cubehelix_palette(len(sorted_atc), start=.5, rot=-.75)
atc_colors = {atc:cor for atc, cor in zip(sorted_atc, _colors)}

# Determine doubling times
stats = stats[((stats['carbon']=='glucose') | (stats['carbon']=='acetate') |
                (stats['carbon']=='glycerol')) & 
                ((stats['temp']==37) |  (stats['temp']==32) | 
                (stats['temp']==42))] 

tidy_stats = pd.DataFrame([])
for g, d in stats.groupby(['date', 'carbon', 'temp', 'run_number']):
    growth_rate = d[d['parameter']=='max df']['value'].values[0]
    growth_err = d[d['parameter']=='max df std']['value'].values[0]
    dbl_time = d[d['parameter']=='inverse max df']['value'].values[0]
    dbl_err = d[d['parameter']=='inverse max df std']['value'].values[0]

    tidy_stats = tidy_stats.append({'date':g[0], 'carbon':g[1], 'temp_C':g[2], 'run_number':g[3],
                                    'growth_rate':growth_rate, 
                                    'dbl_time':dbl_time,
                                    'growth_err':growth_err,
                                    'dbl_err':dbl_err}, 
                                    ignore_index=True)
tidy_stats['growth_rate'] *= 60
tidy_stats['growth_err'] *= 60

# Summarize the growth rates
tidy_stats = tidy_stats.groupby(['carbon', 'temp_C']).agg(('mean', 'sem')).reset_index()
for g, d in tidy_stats.groupby(['carbon', 'temp_C']):
    foldchange.loc[(foldchange['carbon']==g[0]) & (foldchange['temp']==g[1]),
                   'rate_mean'] = d['growth_rate']['mean'].values[0]
    foldchange.loc[(foldchange['carbon']==g[0]) & (foldchange['temp']==g[1]),
                   'rate_sem'] = d['growth_rate']['sem'].values[0]
    foldchange.loc[(foldchange['carbon']==g[0]) & (foldchange['temp']==g[1]),
                   'dbl_mean'] = d['dbl_time']['mean'].values[0]
    foldchange.loc[(foldchange['carbon']==g[0]) & (foldchange['temp']==g[1]),
                   'dbl_sem'] = d['dbl_time']['sem'].values[0]


# Restrict the fold-change measurements to the dilution circuit
fc = foldchange[(foldchange['strain']=='dilution')]
fc = fc[fc['repressors'] >= 10]
#%%

# Set up the kind of complicated figure canvas
fig = plt.figure( figsize=(5, 4.5), dpi=100)
gs = GridSpec(3, 2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1:, 0])
ax3 = fig.add_subplot(gs[0, 1])
ax4 = fig.add_subplot(gs[1:, 1])

bins = np.logspace(0, 4, 75)
high_conc = fc[fc['atc_ngml']==7]
for g, d in high_conc.groupby(['carbon', 'temp']):
    if g[1] == 37:
        ax1.hist(d['repressors'], bins=bins, label=f'{g[0]}, 37°C', density=True,
                edgecolor=edge_colors[g[0]], facecolor=fill_colors[g[0]],
                alpha=0.3)
        ax1.plot(d['repressors'].mean(), 0.005, marker='v', 
                markerfacecolor=fill_colors[g[0]], markeredgecolor=edge_colors[g[0]], alpha=0.75)
    if (g[0] == 'glucose'):
        ax3.hist(d['repressors'], bins=bins, label=f'glucose, {g[1]}°C', density=True,
                edgecolor=edge_colors[g[1]], facecolor=fill_colors[g[1]],
                alpha=0.3)
        print(d['repressors'].mean())

        ax3.plot(d['repressors'].mean(), 0.005, marker='v', 
                markerfacecolor=fill_colors[g[1]], markeredgecolor=edge_colors[g[1]], alpha=0.75)


for g, d in fc.groupby(['atc_ngml']):
    d.sort_values('dbl_mean', inplace=True)
    # Isolate temps
    d_carb = d[d['temp']==37]
    d_temp = d[d['carbon']=='glucose']
    d_carb_grouped = d_carb.groupby('carbon').agg(('mean', 'sem')).reset_index()
    d_carb_grouped.sort_values(('dbl_mean', 'mean'), inplace=True)
    d_temp_grouped = d_temp.groupby('temp').agg(('mean', 'sem')).reset_index()
    d_temp_grouped.sort_values(('dbl_mean', 'mean'), inplace=True)

    ax2.errorbar(d_carb_grouped['dbl_mean']['mean'], d_carb_grouped['repressors']['mean'],
                d_carb_grouped['repressors']['sem'], capsize=1.5, lw=0.5, color=atc_colors[g])
    ax4.errorbar(d_temp_grouped['dbl_mean']['mean'], d_temp_grouped['repressors']['mean'],
                d_temp_grouped['repressors']['sem'], capsize=1.5, lw=0.5, color=atc_colors[g])


for a in [ax1, ax3]:
    a.legend(fontsize=6, handlelength=1)
    a.set_xlim([10, 1500])
    a.set_xscale('log')
    a.set_ylim([0, 0.0065])
    a.set_xlabel('repressors per cell', fontsize=8)
    a.set_ylabel('frequency', fontsize=8)

plt.tight_layout()
 #%%
