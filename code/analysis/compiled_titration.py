#%% 
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import bokeh.io
import bokeh.plotting
import mwc.viz
colors, color_list = mwc.viz.bokeh_theme()
bokeh.io.output_notebook()

#%%
# Load the titration data
data = pd.read_csv('../../data/raw_compiled_snaps.csv')
data.head()

# COnsider only glucose at 37°
CARBON = 'glucose'
TEMP = 37 
STRAIN = 'dilution'
sel = data[(data['carbon']==CARBON) & (data['temp']==TEMP)]
          
# Subtract the autofluorescence for that day.
sub_df = []
for g, d in sel.groupby(['carbon', 'temp', 'date', 'run_number']):
    d = d.copy()

    # Determine the autofluorescence
    auto_yfp = d[d['strain']=='auto']['fluor1_mean_death'].mean()
    auto_mch = d[d['strain']=='delta']['fluor2_mean_death'].mean()
    d = d[d['strain']=='dilution']

    # Perform the subtraction and add to the storage list. 
    d['fluor1_sub'] = d['fluor1_mean_death'] - auto_yfp
    d['fluor2_sub'] = d['fluor2_mean_death'] - auto_mch
    sub_df.append(d)
sub_df = pd.concat(sub_df)

#%%
grouped = sub_df.groupby(['carbon', 'date', 'run_number',
                          'temp', 'atc_ngml']).agg(('mean', 'sem')).reset_index()
                          
# Instantiate the figure
fig, ax = plt.subplots(1, 1, figsize=(2.3, 2.3), dpi=150)
ax.set_title(f'{CARBON}, {TEMP}°C', style='italic', loc='left', fontsize=8)
ax2 = ax.twinx()
ax.set_ylim([0, 150])
ax2.set_ylim([0, 600])

# Add labels
ax.set_xlabel('ATC [ng / mL]', style='italic')
ax.set_ylabel('mean mCherry intensity [a.u. / pixel]', style='italic')
ax2.set_ylabel('mean YFP intensity [a.u. / pixel]', style='italic')
ax2.grid(False)
ax2.spines['right'].set_visible(True)
ax2.spines['right'].set_color('#444147')


for g, d in grouped.groupby(['carbon', 'temp', 'date', 'run_number']):
    if (g[-2] != 20190617) & (g[-2] != 20190611):        
        ax.plot(d['atc_ngml'], d['fluor2_sub']['mean'], '-', 
                lw=0.5, color=colors['dark_red'])             
        ax.errorbar(d['atc_ngml'], d['fluor2_sub']['mean'], d['fluor2_sub']['sem'],
                    fmt='.', ms=5, markerfacecolor=colors['light_red'], 
                    markeredgecolor=colors['dark_red'], capsize=1, lw=0.5,
                    color=colors['dark_red'], markeredgewidth=0.5)

        ax2.plot(d['atc_ngml'], d['fluor1_sub']['mean'], '-', 
                lw=0.5, color=colors['dark_orange'])             
        ax2.errorbar(d['atc_ngml'], d['fluor1_sub']['mean'], d['fluor1_sub']['sem'],
                    fmt='.', ms=5, markerfacecolor=colors['light_orange'], 
                    markeredgecolor=colors['dark_orange'], capsize=1, lw=0.5,
                    color=colors['orange'], markeredgewidth=0.5)

plt.savefig('../figures/example_titration.pdf', bbox_inches='tight', facecolor='white')

#%%


#%%
