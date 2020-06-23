#%%
import numpy as np
import matplotlib.pyplot as plt
import tqdm
import pandas as pd
import mwc.bayes
import mwc.stats
import mwc.process
import scipy.stats
import mwc.viz
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
colors, palette = mwc.viz.plotting_style()

# Load the fluctuation data
data = pd.read_csv('../../data/analyzed_foldchange.csv', comment='#')

# Keep only the dilution strain and the glucose samples as well as >0 reps
data = mwc.process.condition_filter(data, carbon='glucose')
data['repressors'] *= 1.16

# Group by date, run_number, and ATC concentration to compute the mean fc
grouped = data.groupby(['date', 'run_number', 'atc_ngml', 'temp']).mean().reset_index()

# Determine the number of unique temperatures and add an identifier. 
grouped['idx'] = grouped.groupby('temp').ngroup() + 1
summarized = grouped.groupby(['atc_ngml', 'temp'])[['repressors', 'fold_change']].agg(('mean', 'sem'))

#%%
# Load the inferential model. 
model = mwc.bayes.StanModel('../stan/constrained_entropy_estimation.stan', force_compile=True)
epsilon_model = mwc.bayes.StanModel('../stan/DNA_binding_energy.stan', force_compile=True)

#%%
# Assign the data dictionary. 
_samples = []
data_dict = {'Nns': 4.6E6}
for idx in tqdm.tqdm([1, 2, 3]):
    data_dict['J'] =  grouped['idx'].max()
    data_dict['N'] = len(grouped) 
    data_dict['J_REF'] = idx 
    data_dict['idx'] = grouped['idx']
    data_dict['temp'] = np.array([32, 37, 42]) + 273.15
    data_dict['repressors'] = grouped['repressors']
    data_dict['foldchange'] = grouped['fold_change']
    print(data_dict['J_REF'])
    # Sample the model.
    fit, samples = model.sample(data_dict, iter=5000, control=dict(adapt_delta=0.99))

    # Extract and summarize the parameters
    params = model.summarize_parameters()

    # Add identifiers to the samples
    renamed_params = ['epRA_star', 'delta_SR', 'ref_epRA', 'sigma']

    keymap = {}
    for dim, temp in zip(grouped['idx'].unique(), grouped['temp'].unique()):
        params.loc[params['dimension']==dim, 'temp'] = temp
        keymap[dim] = temp
        keymap[temp] = dim

    samples_dfs = []
    for g, d in params.groupby(['temp']):
        for p in renamed_params:
            df = pd.DataFrame()
            df['lp__'] = samples['lp__']
            if (p == 'epRA_star'):
                df['value'] = samples[f'{p}[{keymap[g]}]']   
            else:
                df['value'] = samples[f'{p}']
            df['parameter'] = p
            df['step_id'] = np.arange(len(df)) + 1
            df['temp'] = g
            df['reference_temp'] = data_dict['temp'][idx - 1] - 273.15

            samples_dfs.append(df)
    samples = pd.concat(samples_dfs)
    _samples.append(samples)

lf_samples = pd.concat(_samples)
lf_samples.to_csv('../../data/constrained_entropy_estimation_samples.csv', index=False)

# Infer te DNA binding energies
#%%
dfs = []
for g, d in tqdm.tqdm(grouped.groupby(['temp'])):
    data_dict = {'N':len(d), 'repressors':d['repressors'],
                'Nns':4.6E6, 'foldchange':d['fold_change']}
    fit, samples = epsilon_model.sample(data_dict, iter=5000, 
                                        control=dict(adapt_delta=0.99))
    params = epsilon_model.summarize_parameters()
    params['temp'] = g
    dfs.append(params)

epsilon_params = pd.concat(dfs, sort=False)

epsilon_params.to_csv('../../data/temperature_DNA_binding_energy_summary.csv',
                        index=False)
#%%
# Define constants for the curves
rep_range = np.logspace(0, 3, 200)
temp_range = np.linspace(30, 45, 200) + 273.15

temp_colors = {32:colors['blue'], 37: colors['purple'], 42:colors['red']}

# Set up the figure canvases
fig, ax = plt.subplots(2, 3, figsize=(6, 4))
for a in ax:
    mwc.viz.despine(a)

idx_map = {1:32, 2:37, 3:42}
for i in range(3):
    ax[0, i].set_xscale('log')
    ax[0, i].set_xlabel('repressors per cell')
    ax[0, i].set_yscale('log')
    ax[0, i].set_ylabel('fold-change')
    ax[1, i].set_xlabel('1000 / T [K$^{-1}$]')
    ax[1, i].set_ylabel(r'$\varepsilon$ [$k_BT$]')
    for g, d in epsilon_params.groupby(['temp']):
        if g == idx_map[i + 1]:
            face = colors['grey']
            edge = temp_colors[g]
            ms = 3.5
        else:
            face = temp_colors[g]
            edge = colors['grey']
            ms=4
        ax[1, i].vlines(1000/(g + 273.15), d[d['parameter']=='epRA']['hpd_min'],
                        d[d['parameter']=='epRA']['hpd_max'], lw=0.5,
                        color=temp_colors[g], label='__nolegend__')

        ax[1, i].plot(1000/(g + 273.15), d[d['parameter']=='epRA']['median'], 'o',
                    color=temp_colors[g], label=g, ms=ms, markerfacecolor=face,
                    markeredgecolor=edge, markeredgewidth=0.5)

# Generate the van't hoff plot        
temp_range = np.linspace(30, 45, 200)
for i, t in enumerate(idx_map.values()):
    _ref = lf_samples[lf_samples['reference_temp']==t]
    delS = _ref[_ref['parameter']=='delta_SR']['value'].unique()
    ref_epRA = _ref[(_ref['parameter']=='epRA_star') & (_ref['temp']==t)]['value'].unique()
    vant_hoff_cred = np.zeros((2, len(temp_range))) 
    medians = []
    for j, _t in enumerate(temp_range):
        _vant = delS * (t - _t) + ref_epRA
        vant_hoff_cred[:, j] = mwc.stats.compute_hpd(_vant, 0.95)
        medians.append(np.median(_vant))

    ax[1, i].fill_between(1000/(temp_range + 273.15), vant_hoff_cred[0, :], 
                        vant_hoff_cred[1, :], color=temp_colors[t], alpha=0.25,
                        label='__nolegend__')
    ax[1, 0].plot(1000/(temp_range + 273.15), medians, '--', color=temp_colors[t], lw=0.5, 
                        label='__nolegend__', alpha=0.75)

    for g, d in summarized.groupby(['temp']):
        epra = _ref[(_ref['temp']==g) & (_ref['parameter']=='epRA_star')]['value'].values
        if g == t:
            face = colors['grey']
            edge = temp_colors[g] 
            ms=3.5
        else:
            face = temp_colors[g]
            edge = colors['grey']
            ms=4
        ax[0, i].errorbar(d['repressors']['mean'], d['fold_change']['mean'],
                          xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'],
                          fmt='o', linestyle='none', lw=0.5, capsize=1,
                          markerfacecolor=face, markeredgecolor=edge, 
                          markeredgewidth=0.5, label=g, color=temp_colors[g], 
                          ms=ms)
        cred_region = np.zeros((2, len(rep_range)))

        for j, r in enumerate(rep_range):
            arch = (1 + (r/4.6E6) * np.exp(-epra))**-1
            cred_region[:,  j] = mwc.stats.compute_hpd(arch, 0.95)
        ax[0, i].fill_between(rep_range, cred_region[0, :], cred_region[1, :],
                            color=temp_colors[g], alpha=0.45, label='__nolegend__')
plt.tight_layout()
# %%
# Define the supplemental figure
fig, ax = plt.subplots(1, 3, figsize=(6, 2))
mwc.viz.despine(ax)
ax[0].set_xlabel('repressors per cell')
ax[0].set_ylabel('fold-change')
ax[1].set_xlabel(r'$\frac{1000}{T}$ [K$^{-1}$]')
ax[1].set_ylabel('effective DNA binding\nenergy ε [$k_BT$]')
ax[1].set_ylim([-17, -12])
ax[2].set_xlabel(r'$\Delta S$ [$k_B$]')
ax[2].set_ylabel(r'cumulative distribution')
ax[0].set_xscale('log')
ax[0].set_yscale('log')

# Isolate the parameters
_ref = lf_samples[lf_samples['reference_temp']==37]
ref_epra = _ref[_ref['parameter']=='ref_epRA']['value'].unique()
delS = _ref[_ref['parameter']=='delta_SR']['value'].unique()

# Define the repressor range
rep_range = np.logspace(0, 3, 200)

for g, d in summarized.groupby(['temp']):

    # plot the fold-chagne curves 
    cred_region = np.zeros((2, len(rep_range)))
    for i, r in enumerate(rep_range):
        fc = (1 + (r/4.6E6) * np.exp(-(delS * (37 - g) + ref_epra)))**-1
        cred_region[:, i] = mwc.stats.compute_hpd(fc, 0.95)
    ax[0].fill_between(rep_range, cred_region[0, :], cred_region[1, :], color=temp_colors[g],
                    label='__nolegend__', alpha=0.5)
    if g == 37:
        face = colors['grey']
        edge = temp_colors[g]
        ms = 3.5
    else:
        face = temp_colors[g]
        edge = colors['grey']
        ms = 4
    ax[0].errorbar(d['repressors']['mean'], d['fold_change']['mean'],
                   xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'], 
                   ms=ms, color=temp_colors[g], lw=0.5, markeredgewidth=0.5,
                   markerfacecolor=face, markeredgecolor=edge, label=f'{g}° C',
                   fmt='o', linestyle='none')

# Plot the van't hoff curve. 
temp_range = np.linspace(30, 45, 100)
vant_cred = np.zeros((2, len(temp_range)))

for i, t in enumerate(temp_range):
    vant_cred[:, i] = mwc.stats.compute_hpd(delS * (37 - t) + ref_epra, 0.95)

ax[1].fill_between(1000/(temp_range + 273.15), vant_cred[0, :], vant_cred[1, :],
                   color=colors['light_grey'], alpha=0.25, label=r'$-\frac{\Delta H}{k_BT} + \frac{\Delta S}{k_B}$')

# Generate the van't hoff plot
for g, d in epsilon_params.groupby(['temp']):
    face = temp_colors[g]
    edge = colors['grey']
    ms=4
    ax[1].vlines(1000/(g + 273.15), d[d['parameter']=='epRA']['hpd_min'],
                    d[d['parameter']=='epRA']['hpd_max'], lw=0.5,
                    color=temp_colors[g], label='__nolegend__')

    ax[1].plot(1000/(g + 273.15), d[d['parameter']=='epRA']['median'], 'o',
                    color=temp_colors[g], label=f'{g}° C', ms=ms, markerfacecolor=face,
                    markeredgecolor=edge, markeredgewidth=0.5)
for t in [32, 37, 42]:
    _ref = lf_samples[lf_samples['reference_temp']==t]
    delS = _ref[_ref['parameter']=='delta_SR']['value'].unique()
    x, y = np.sort(delS), np.arange(0, len(delS), 1) / len(delS)
    ax[2].step(x, y, lw=0.75, alpha=0.75, color=temp_colors[t], label=f'{t}° C')
plt.tight_layout()
for a in ax:
    a.legend(fontsize=5.5)

fig.text(0.01, 0.95, '(A)', fontsize=8)
fig.text(0.33, 0.95, '(B)', fontsize=8)
fig.text(0.69, 0.95, '(C)', fontsize=8)
plt.savefig('../../figs/vant_hoff_analysis.pdf')
# %%
