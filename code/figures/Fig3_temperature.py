#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import mwc.model
import mwc.process
import mwc.stats
colors, palette = mwc.viz.plotting_style()
constants = mwc.model.load_constants()

# Load the experimental measurements
data = pd.read_csv('../../data/analyzed_foldchange.csv', comment='#')
data = mwc.process.condition_filter(data, carbon='glucose')
data['repressors'] *= 1.16

# Load the data collapse data
old = pd.read_csv('../../data/Garcia2011_Brewster2014_RazoMejia2018_Chure2019.csv', 
                  comment='#')

# Load the sampling information. 
entropy = pd.read_csv('../../data/pooled_entropic_parameter_samples.csv')

#%%
# Set up the figure canvas
fig, ax = plt.subplots(2, 2, figsize=(5, 5))

# Format the axes and add titles
mwc.viz.despine(ax.ravel())
for a in ax.ravel():
    a.set_yscale('log')
    a.set_ylabel('fold-change')
titles = ['rescaling of thermal energy', 'inclusion of entropic penalty']

for i in range(2):
    ax[0, i].set_ylim([0.005, 1.1])
    ax[0, i].set_xlim([1, 1E3])
    ax[1, i].set_xlim([-5, 5])
    ax[1, i].set_ylim([5E-3, 1.5])
    ax[0, i].set_xscale('log')
    ax[0, i].set_xlabel('repressors per cell')
    ax[1, i].set_xlabel('free energy [$k_BT$]')
    mwc.viz.titlebox(ax[0, i], titles[i], color=colors['black'])

# Plot the master curve and collapse function
bohr_range = np.linspace(-10, 10, 200)
master_curve = (1 + np.exp(-bohr_range))**-1
for i in range(2):
    ax[1, i].plot(old['bohr_parameter'], old['mean'], '.', label='previous data',
                ms=5, color='lightgrey')
    ax[1, i].plot(bohr_range, master_curve, 'k-', lw=0.75, label='scaling function')



# Compute the summarized data
summarized = data.groupby(['atc_ngml', 'date', 'temp']).mean().reset_index()
summarized = summarized.groupby(['atc_ngml', 
                                'temp'])[
                               ['repressors', 'fold_change']
                               ].agg(('mean', 'sem')).reset_index()
# Define the colors
temp_colors = {37:colors['purple'], 32:colors['blue'], 42:colors['red']}
for g, d in summarized.groupby(['temp']):
   for i in range(2):
       ax[0, i].errorbar(d['repressors']['mean'], d['fold_change']['mean'],
                        xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'],
                        fmt='o', ms=3.75, color=temp_colors[g],
                        markeredgewidth=0.5, markeredgecolor=colors['grey'],
                        lw=0.75, linestyle='none', label=f'glucose, {g} °C') 


# Plot the 37C predictions
ep_r = [-14.1, -13.7]
rep_range = np.logspace(0, 3, 200)
ep, r = np.meshgrid(ep_r, rep_range)
fc = mwc.model.SimpleRepression(r, ep, ka=constants['ka'], ki=constants['ki'],
                                ep_ai=constants['ep_ai'], effector_conc=0).fold_change()
for i in range(2):
    ax[0, i].fill_between(rep_range, fc[:, 0], fc[:, 1], color=temp_colors[37],
                          alpha=0.3)

# Plot the 37C data collapse. 
data_37 = summarized[summarized['temp']==37]
bohr = -mwc.model.SimpleRepression(data_37['repressors']['mean'], ep_r=-13.9,
                                 ka=constants['ka'], ki=constants['ki'], 
                                 ep_ai=constants['ep_ai'], effector_conc=0).bohr_parameter()
for i in range(2):
    ax[1, i].errorbar(bohr, data_37['fold_change']['mean'],
                      data_37['fold_change']['sem'], color=temp_colors[37],
                      markeredgecolor=colors['grey'], markeredgewidth=0.5,
                      fmt='o', linestyle='none', lw=0.75, ms=3.5, label=f'glucose, 37°C',
                      zorder=1000)

# Plot the predictions for simple reweighting. 
T_ref = 37 + 273.15
for t in [32, 42]:
    rel_T = T_ref / (t + 273.15) 
    r, ep = np.meshgrid(rep_range, np.array(ep_r) * rel_T)
    fc = mwc.model.SimpleRepression(r, ep, ka=constants['ka'], 
                                   ki=constants['ki'], ep_ai=rel_T * constants['ep_ai'],
                                   effector_conc=0).fold_change()
    ax[0, 0].fill_between(rep_range, fc[0, :], fc[1, : ], color=temp_colors[t],
                          alpha=0.3, label='__nolegend__')

    # Plot the bohr parameters for the simple reweighting. 
    _d = summarized[summarized['temp']==t]
    bohr = -mwc.model.SimpleRepression(_d['repressors']['mean'], ep_r=rel_T * -13.9,
                                      ka=constants['ka'], ki=constants['ki'],
                                      ep_ai=constants['ep_ai'] * rel_T,
                                      effector_conc=0).bohr_parameter()
    ax[1, 0].errorbar(bohr, _d['fold_change']['mean'], _d['fold_change']['sem'],
                    fmt='o', ms=3.5, markeredgewidth=0.5, markeredgecolor=colors['grey'],
                    color=temp_colors[t], label=f'glucose, {t}°C', zorder=1000)


# Plot the predictions for the inclusion of the entropic parameter. 
for t in [32, 42]:
    cred_region = np.zeros((2, len(rep_range)))
    _sample = entropy[entropy['temp']==t]
    epRA_star = _sample[_sample['parameter']=='epRA_star']['value'].values
    epAI_star = _sample[_sample['parameter']=='epAI_star']['value'].values
    for i, r in enumerate(rep_range):
        fc = mwc.model.SimpleRepression(r, epRA_star, ka=constants['ka'],   
                                        ki=constants['ki'], ep_ai=epAI_star,
                                        effector_conc=0).fold_change() 
        cred_region[:, i] = mwc.stats.compute_hpd(fc, 0.95) 
    ax[0, 1].fill_between(rep_range, cred_region[0, :], cred_region[1, :],  
                        color=temp_colors[t], alpha=0.3, label='__nolegend__')

    # Plot the bohr parameter
    _d = summarized[summarized['temp']==t]
    bohr = -mwc.model.SimpleRepression(R=_d['repressors']['mean'], ep_r=np.median(epRA_star),
                                      ka=constants['ka'], ki=constants['ki'],
                                      ep_ai=np.median(epAI_star), effector_conc=0).bohr_parameter()
    ax[1, 1].errorbar(bohr, _d['fold_change']['mean'],  _d['fold_change']['sem'],
                      fmt='o', ms=3.5, color=temp_colors[t], markeredgewidth=0.5,
                      markeredgecolor=colors['grey'], lw=0.75, linestyle='none',
                      label=f'glucose, {t}°C', zorder=1000)

for a in ax.ravel():
    a.legend(fontsize=6)
fig.text(0.015, 0.99, '(A)', fontsize=8)
fig.text(0.515, 0.99, '(B)', fontsize=8)
plt.tight_layout()
plt.savefig('../../figs/Fig3_temperature.pdf', bbox_inches='tight')
# %%
