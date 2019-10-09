# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import mwc.viz
import seaborn as sns
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
_colors = sns.color_palette('viridis', n_colors=len(sorted_atc)) # + 1)
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
fc = fc[fc['repressors'] >=10]

# Set up the kind of complicated figure canvas
fig = plt.figure( figsize=(5, 5), dpi=150)
gs = GridSpec(7, 4)

# Add the ATC induction plots
ax1 = fig.add_subplot(gs[0:3, 0:2])
ax2 = fig.add_subplot(gs[4:, 0:2])
ax3 = fig.add_subplot(gs[0:3, 2:4])
ax4 = fig.add_subplot(gs[4:, 2:4])



i = 1
for g, d in fc.groupby(['atc_ngml']):
    d.sort_values('dbl_mean', inplace=True)
    if i%1 == 0:
        # Isolate temps
        d_carb = d[d['temp']==37]
        d_temp = d[d['carbon']=='glucose']
        d_carb_grouped = d_carb.groupby(['carbon', 'date', 'run_number']).mean()
        d_carb_grouped = d_carb_grouped.groupby(['carbon']).agg(('mean', 'sem'))
        d_carb_grouped.sort_values(('rate_mean', 'mean'), inplace=True)
        d_temp_grouped = d_temp.groupby(['temp', 'date', 'run_number']).mean()
        d_temp_grouped = d_temp_grouped.groupby('temp').agg(('mean', 'sem')).reset_index()
        d_temp_grouped.sort_values(('rate_mean', 'mean'), inplace=True)

        ax2.errorbar(d_carb_grouped['rate_mean']['mean'], d_carb_grouped['repressors']['mean'],
                    d_carb_grouped['repressors']['sem'], capsize=2, lw=0.5, color=atc_colors[g],
                    fmt='.', linestyle='-', label=g,
                    markeredgecolor='k', markeredgewidth=0.25)
        ax4.errorbar(d_temp_grouped['rate_mean']['mean'], d_temp_grouped['repressors']['mean'],
                    d_temp_grouped['repressors']['sem'], capsize=2, lw=0.5, color=atc_colors[g],
                    fmt='.', linestyle='-', markeredgecolor='k', markeredgewidth=0.25)
        ax1.errorbar(d_carb_grouped['rate_mean']['mean'], d_carb_grouped['area']['mean'],
                    d_carb_grouped['area']['sem'], capsize=2, lw=0.5, color=atc_colors[g],
                    fmt='.', linestyle='-', label=g,
                    markeredgecolor='k', markeredgewidth=0.25)
        ax3.errorbar(d_temp_grouped['rate_mean']['mean'], d_temp_grouped['area']['mean'],
                    d_temp_grouped['area']['sem'], capsize=2, lw=0.5, color=atc_colors[g],
                    fmt='.', linestyle='-', markeredgecolor='k', markeredgewidth=0.25)
 
    i += 1

for a in [ax1, ax3]:
    a.set_ylim([1.6, 2.8])
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8, style='italic')
    a.set_ylabel('projected cell area [µm$^2$]', fontsize=8, style='italic')

for a in [ax2, ax4]:
    a.set_ylabel('repressors per cell', fontsize=8, style='italic')
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8, style='italic')

for a in [ax1, ax2]:
    a.set_xlim([0.18, 0.65])

for a in [ax3, ax4]:
    a.set_xlim([0.44, 0.65])

handles, labels = ax2.get_legend_handles_labels()
leg = ax3.legend(reversed(handles), reversed(labels), title='   ATC\n[ng / mL]', 
                fontsize=6, bbox_to_anchor=(1.15, 1))
leg.get_title().set_fontsize(6)
plt.subplots_adjust(hspace=0.5, wspace=0.6)
plt.savefig('../../figs/Fig_expression_scaling.svg', bbox_inches='tight', 
            facecolor='white')

%il%
