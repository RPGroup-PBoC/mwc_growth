#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import mwc.model
colors, _ = mwc.viz.personal_style()


#%% Load the data sets
ind_data = pd.read_csv('../../data/RazoMejia_2018.csv', comment='#')
old_gods = pd.read_csv('../../data/Garcia2011_Brewster2014.csv', comment='#')
data = pd.read_csv('../../data/analyzed_foldchange.csv', comment='#')

# Restrict to only O2 and prune
ind_data = ind_data[(ind_data['repressors'] > 0) & (ind_data['operator']=='O2') &
                    (ind_data['IPTG_uM']==0)]
ind_data['repressors'] *= 2
ind_data.rename(columns={'fold_change_A':'fold_change'}, inplace=True)
old_gods = old_gods[old_gods['operator']=='O2']
data = data[(data['strain']=='dilution') & (data['carbon']=='glucose') & 
            (data['temp']==37) & (data['repressors'] > 0) & (data['fold_change'] >= 0)]

# Group and summarize as needed
ind_data = ind_data.groupby('repressors').agg(('mean', 'sem')).reset_index()
data = data.groupby(['date', 'run_number', 'atc_ngml']).mean().reset_index()
data = data.groupby('atc_ngml').agg(('mean', 'sem'))

# %%
# Define constants for plotting
rep_range = np.logspace(0, 4, 200)
epR = -13.9
epAI = 4.5
theo = mwc.model.SimpleRepression(R=rep_range, ep_r=epR, ka=139, ki=0.53, ep_ai=epAI,
                                effector_conc=0).fold_change()


# %%
# instantiate the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(3, 3))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('repressors per cell')
ax.set_ylabel('fold-change')
ax.set_xlim([1, 2000])

# Plot the theory.
ax.plot(rep_range, theo, 'k-', zorder=100)

# Plot the induction measurements
ax.errorbar(ind_data['repressors'], ind_data['fold_change']['mean'], 
            ind_data['fold_change']['sem'], fmt='.', ms=4, capsize=1, lw=0.75,
            color=colors['red'], label='Razo-Mejia et al. 2018')

# Plot the old gods measurements. 
garcia = old_gods[old_gods['author']=='garcia']
brewster = old_gods[old_gods['author']=='brewster']
ax.plot(garcia['repressor'], garcia['fold_change'], '.', ms=4, color=colors['blue'], 
        label='Garcia & Phillips, 2011')
ax.plot(brewster['repressor'], brewster['fold_change'], '.', ms=4, color=colors['brown'], 
         label='Brewster et al. 2014')

# Plot the measurements with the correction factor
ax.errorbar(data['repressors']['mean']/2, data['fold_change']['mean'], 
            xerr=data['repressors']['sem']/2, yerr=data['fold_change']['sem'],
            fmt='.', color=colors['light_purple'], markerfacecolor=colors['grey'],
            capsize=1, linestyle='none', label='This work (without correction)',
            zorder=101, ms=6, alpha=0.75)

ax.errorbar(data['repressors']['mean'], data['fold_change']['mean'], 
            xerr=data['repressors']['sem'], yerr=data['fold_change']['sem'],
            fmt='.', color=colors['dark_purple'], markerfacecolor=colors['light_purple'],
            capsize=1, linestyle='none', label='This work (with correction)',
            zorder=102, ms=6)

ax.legend()
plt.savefig('../../figs/FigSX_cal_factor_correction.pdf', bbox_inches='tight', 
            facecolor='white')
# %%
