#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import mwc.model
import mwc.process
colors, palette = mwc.viz.plotting_style()
constants = mwc.model.load_constants()

# Load the data sets. 
data = pd.read_csv('../../data/analyzed_foldchange.csv', comment='#')
data['repressors'] *= 1.16
carbon_data = mwc.process.condition_filter(data, temp=37)
temp_data = mwc.process.condition_filter(data, carbon='glucose')
entropy = pd.read_csv('../../data/pooled_entropic_parameter_samples.csv', comment='#')
old = pd.read_csv('../../data/Garcia2011_Brewster2014_RazoMejia2018_Chure2019.csv', comment='#')
old = old[old['repressors'] > 0 ]

#%% Set up the figure canvas
fig, ax = plt.subplots(2, 4, figsize=(7, 4), sharey=True)
mwc.viz.despine(ax.ravel())

# Format the axes as necessary
ax[0, 0].set_ylabel('fold-change')
ax[1, 0].set_ylabel('fold-change')
ax[0, 1].set_xlabel('IPTG [µM]')

for a in [ax[0, 0], ax[0, 2], ax[0, 3]]:
    a.set_xscale('log')
    a.set_xlabel('repressors per cell') 
    a.set_yscale('log')

for i in range(4):
    ax[1, i].set_xlabel('free energy [$k_BT$]')
    ax[1, i].set_yscale('log')
    ax[0, i].set_xscale('log')

# Plot the master curve for the collapse 
bohr = np.linspace(-10, 10, 200)
collapse = (1 + np.exp(-bohr))**-1
for i in range(4):
    ax[1, i].plot(bohr, collapse, 'k-', label='scaling function')

# Plot the garcia and brewster data. 
auth_glyphs = {'garcia': 'o', 'brewster': 's', 'razo-mejia': 'D'}
op_colors = {'Oid':colors['light_grey'], 'O1':colors['blue'], 
              'O2':colors['purple'], 'O3':colors['green']}
for g, d in old[old['IPTGuM']==0].groupby(['author', 'operator', 'repressors']):
    if g[0] == 'chure':
        continue
    if (np.round(g[-1], decimals=-1) == 260) & (g[0]=='garcia'):
        face = 'white'
        edge = op_colors[g[1]] 
        zorder = 1000
    else:
        face = op_colors[g[1]]
        edge = colors['grey']
        zorder = 999
    
    ax[0, 0].plot(d['repressors'].mean(), d['mean'].mean(), linestyle='none', marker=auth_glyphs[g[0]],
                  color=face, markeredgecolor=edge, markeredgewidth=0.5, ms=3.5,
                  zorder=zorder)
    ax[1, 0].plot(d['bohr_parameter'].mean(), d['mean'].mean(), linestyle='none', marker=auth_glyphs[g[0]],
                color=face, markeredgecolor=edge, markeredgewidth=0.5, ms=3.5,
                zorder=zorder)

# Plot the leakiness predictions
rep_range = np.logspace(0, 3.5, 200)
for o in ['O1', 'O2', 'O3', 'Oid']:
    fc = mwc.model.SimpleRepression(R=rep_range, ep_r=constants[o],
                                    ka=constants['ka'], ki=constants['ki'],
                                    ep_ai=constants['ep_ai'], effector_conc=0).fold_change()
    ax[0, 0].plot(rep_range, fc, color=op_colors[o], lw=0.75)


# Plot the induction data. 
op_glyphs = {'O1':'s', 'O2':'o', 'O3':'D'}
rep_colors = {22:colors['blue'], 60:colors['green'], 124:colors['purple'],
              260:colors['orange'], 1220:colors['red'], 1740:colors['brown']}

for g, d in old[old['author']=='razo-mejia'].groupby(['operator', 'repressors', 'IPTGuM']):
    if (g[0] == 'O2') & (g[1] == 260):
        edge = rep_colors[g[1]]
        face = 'white' 
        zorder = 1000
    else:
        edge = colors['grey']
        face = rep_colors[g[1]]
        zorder = 999

    ax[0, 1].plot(d['IPTGuM'].mean(), d['mean'].mean(), linestyle='none', 
                marker=op_glyphs[g[0]], color=face, markeredgewidth=0.5,
                markeredgecolor=edge, ms=3.5, zorder=zorder)
    ax[1,1].plot(d['bohr_parameter'].mean(), d['mean'].mean(), linestyle='none',
                marker=op_glyphs[g[0]], color=face, markeredgewidth=0.5,
                markeredgecolor=edge, ms=3.5, zorder=zorder)
            
# Plot the induction profiles 
c_range = np.logspace(-2, 4, 300)
for o in ['O1', 'O2', 'O3']:
    for r, c in rep_colors.items():
        fc = mwc.model.SimpleRepression(R=r, ep_r=constants[o], ka=constants['ka'],
                                        ki=constants['ki'], ep_ai=constants['ep_ai'],
                                        effector_conc=c_range).fold_change()
        ax[0, 1].plot(c_range, fc, '-', lw=0.75, color=c, label='__nolegend__')

# Plot the carbon source data
carbon_summary = carbon_data.groupby(['atc_ngml', 'carbon', 'date']).mean().reset_index()
carbon = carbon_summary.groupby(['carbon', 'atc_ngml']).mean().reset_index()
carbon_colors = {'glucose':colors['purple'], 
                'glycerol':colors['green'],
                'acetate':colors['brown']}
for g, d in carbon.groupby(['carbon']):
    ax[0, 2].plot(d['repressors'], d['fold_change'], 'o', 
                 markerfacecolor=carbon_colors[g], markeredgecolor=colors['grey'],
                 markeredgewidth=0.5, label='__nolegend__', ms=3.5)

    # Compute the bohr parameter. 
    bohr = -mwc.model.SimpleRepression(R=d['repressors'], ep_r=constants['O2'], 
                                      ka=constants['ka'], ki=constants['ki'],
                                      ep_ai=constants['ep_ai'], effector_conc=0).bohr_parameter()
    ax[1, 2].plot(bohr, d['fold_change'], 'o', ms=3.5, color=carbon_colors[g], 
                markeredgecolor=colors['grey'], markeredgewidth=0.5, zorder=1000)

# Plot the carbon prediction curve. 
fc = mwc.model.SimpleRepression(R=rep_range, ep_r=constants['O2'], 
                                ka=constants['ka'], ki=constants['ki'],
                                ep_ai=constants['ep_ai'], effector_conc=0).fold_change()
ax[0, 2].plot(rep_range, fc, color=colors['purple'], lw=0.75, label='__nolegend__')
                               
# Plot the temperature data 
temp_summary = temp_data.groupby(['atc_ngml', 'date', 'temp']).mean()
temp = temp_summary.groupby(['atc_ngml', 'temp']).mean()
temp_colors = {37:colors['purple'], 32:colors['blue'], 42:colors['red']}
for g, d in temp.groupby(['temp']):
    ax[0, 3].plot(d['repressors'], d['fold_change'], 'o', ms=3.5, color=temp_colors[g],
                  markeredgecolor=colors['grey'], markeredgewidth=0.5, label='__nolegend__')

    # Load the appropriate statistics and compute the predictions
    if g == 37:
        epRA = constants['O2']
        epAI = constants['ep_ai']
    else:
        _summary = entropy[entropy['temp']==g]
        epRA = _summary[_summary['parameter']=='epRA_star']['value'].median()
        epAI = _summary[_summary['parameter']=='epAI_star']['value'].median()
    fc = mwc.model.SimpleRepression(rep_range, ep_r=epRA, ka=constants['ka'],
                                      ki=constants['ki'], ep_ai=epAI, effector_conc=0).fold_change()
    ax[0, 3].plot(rep_range, fc, '-', color=temp_colors[g], lw=0.75, label='__nolegend__')

    # Compute the bohr parameter
    bohr = -mwc.model.SimpleRepression(d['repressors'], epRA, ka=constants['ka'], 
                                      ki=constants['ki'], ep_ai=epAI, effector_conc=0).bohr_parameter()
    ax[1, 3].plot(bohr, d['fold_change'], 'o', ms=3.5, color=temp_colors[g],
                markeredgecolor=colors['grey'], markeredgewidth=0.5, label='__nolegend__')

plt.subplots_adjust(wspace=0.15, hspace=0.35)
plt.savefig('../../figs/Fig4_collapse_contexts.svg', bbox_inches='tight')
# %%


# %%
