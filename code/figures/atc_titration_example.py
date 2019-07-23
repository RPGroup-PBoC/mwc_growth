#%%
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import phd.viz
import mwc.stats
import mwc.viz
color = phd.viz.phd_style()

#%%
# Load the lineages and isolate the glucose samples
snaps = pd.read_csv('../../data/raw_compiled_snaps.csv')
glucose = snaps[(snaps['carbon']=='glucose') & (snaps['temp']==37)]

dfs = []

for g, d in glucose.groupby(['date', 'run_number']):
    d = d.copy()
    # Compute the median auto fluorescence 
    auto_yfp = d[d['strain']=='auto']['fluor1_mean_death'].median()
    auto_mch = d[d['strain']=='auto']['fluor2_mean_death'].median()

    # Subtract autogluorescence and add the dilution strain to the df
    d['yfp_sub'] = (d['fluor1_mean_death'] - auto_yfp) * d['area_death']
    d['mch_sub'] = (d['fluor2_mean_death'] - auto_mch) * d['area_death']
    dfs.append(d[d['strain']=='dilution'])
dilution = pd.concat(dfs)

# %%
# Group by the atc concentration
summary = dilution.groupby(['atc_ngml']).agg(('mean', 'sem')).reset_index()
summary['yfp'] = summary['yfp_sub']['mean'] / summary['yfp_sub']['mean'].max()
summary['yfp_err'] = summary['yfp_sub']['sem'] / summary['yfp_sub']['mean'].max()
summary['mch'] = summary['mch_sub']['mean'] / summary['mch_sub']['mean'].max()
summary['mch_err'] = summary['mch_sub']['sem'] / summary['mch_sub']['mean'].max()

# Set up the figure. 
fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5))
ax.xaxis.set_tick_params(labelsize=7)
ax.yaxis.set_tick_params(labelsize=7)
ax.set_xlabel('ATC [ng / mL]', style='italic', fontsize=8)
ax.set_ylabel('relative fluorescence', style='italic', fontsize=8)
ax.errorbar(summary['atc_ngml'], summary['yfp'], summary['yfp_err'], linestyle='-', fmt='.',
            color=colors['orange'], label='YFP', markerfacecolor=colors['light_orange'],
            capsize=1.5, lw=0.75, ms=6, markeredgewidth=0.5)
ax.errorbar(summary['atc_ngml'], summary['mch'], summary['mch_err'], linestyle='-', fmt='.',
            color=colors['dark_red'], label='mCherry', markerfacecolor=colors['light_red'],
            capsize=1.5, lw=0.75, ms=6, markeredgewidth=0.5)
ax.set_xscale('log')
ax.set_title('glucose, 37 °C', style='italic', loc='left', fontsize=8)
_ = ax.legend(fontsize=8)
plt.savefig('../../figs/Fig2_relative_fluorescence.pdf', bbox_inches='tight',
            facecolor='none')

#%%


#%%
