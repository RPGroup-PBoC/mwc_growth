
#%%
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import phd.viz
import mwc.stats
import mwc.viz
colors = phd.viz.phd_style()

#%%
# Load the lineages and isolate the glucose samples
carbon = 'glucose'
temp = 37 
snaps = pd.read_csv('../../data/raw_compiled_snaps.csv')
# snaps.loc[snaps['atc_ngml']<1, 'atc_ngml'] *= 10
glucose = snaps[(snaps['carbon']==carbon) & (snaps['temp']==temp)]

dfs = []

for g, d in glucose.groupby(['date', 'run_number']):
    d = d.copy()
    # Compute the median auto fluorescence 
    auto_yfp = d[d['strain']=='auto']['fluor1_mean_death'].mean()
    auto_mch = d[d['strain']=='auto']['fluor2_mean_death'].mean()

    # Subtract autogluorescence and add the dilution strain to the df
    d['yfp_sub'] = (d['fluor1_mean_death'] - auto_yfp)
    d['mch_sub'] = (d['fluor2_mean_death'] - auto_mch) 
    dfs.append(d[d['strain']=='dilution'])
dilution = pd.concat(dfs)

# %%
# Group by the atc concentration
n = dilution.groupby(['date', 'run_number']).ngroups

# Instantiate figure
fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5), sharex=True, sharey=True)

ax2 = ax.twinx()
ax2.yaxis.grid(False)
ax2.spines['right'].set_visible(True)
ax.set_xscale('log')
# ax.set_xlim([.1, 11])
ax.xaxis.set_tick_params(labelsize=7)
ax.yaxis.set_tick_params(labelsize=7)
ax.set_xlabel('ATC [ng / mL]', style='italic', fontsize=8)
ax.set_ylabel('YFP intensity [a.u. / pixel]', style='italic', fontsize=8)
ax2.set_ylabel('mCherry intensity [a.u. / pixel]', style='italic', fontsize=8, rotation=-90,
               labelpad=10)
    


summary = dilution.groupby(['atc_ngml', 'date', 'run_number']).agg(('mean', 'sem')).reset_index()
i = 0
# ax = ax.ravel()
x = []
m = []
y = []
for g, d in summary.groupby(['atc_ngml']):
            d = d.copy() 
        # if i == 6:
            print(g)
            d['yfp'] = d['yfp_sub']['mean'] #/ d['yfp_sub']['mean'].min()
            # d.loc[d['atc_ngml'] < 1, 'yfp'] = d.loc[d['atc_ngml'] < 1, 'yfp'].values * 10
            d['yfp_err'] = d['yfp_sub']['sem'] #/ d['yfp_sub']['mean'].min()
            # d.loc[d['atc_ngml'] < 1, 'yfp_err']  *= 10
            d['mch'] = d['mch_sub']['mean'] #/ d['mch_sub']['mean'].min()
            # d.loc[d['atc_ngml'] < 1, 'mch']  *= 10
            d['mch_err'] = d['mch_sub']['sem'] #/ d['mch_sub']['mean'].min()
            # d.loc[d['atc_ngml'] < 1, 'mch_err']  *= 10
            x.append(g)
            m.append(d['mch'].mean())
            y.append(d['yfp'].mean())
            ax.errorbar(d['atc_ngml'].mean(), d['yfp'].mean(), d['yfp'].std() / np.sqrt(len(d)), linestyle='-', fmt='.',
                    color=colors['orange'], label='__nolegend__', markerfacecolor=colors['light_orange'],
                    capsize=1.5, lw=0.75, ms=6, markeredgewidth=0.5)
            ax2.errorbar(d['atc_ngml'].mean(), d['mch'].mean(), d['mch'].std() / np.sqrt(len(d)), linestyle='-', fmt='.',
                    color=colors['dark_red'], label='__nolegend__', markerfacecolor=colors['light_red'],
                    capsize=1.5, lw=0.75, ms=6, markeredgewidth=0.5)
ax.set_title(f'{carbon}, {temp}° C', style='italic', loc='left', fontsize=8)
        # i+=1
ax2.plot(x, m, color=colors['red'])
ax.plot(x, y, color=colors['orange'])
plt.tight_layout()
plt.savefig(f'../../figs/Fig2_atc_titration.pdf', bbox_inches='tight',
            facecolor='none')

#%%


#%%
