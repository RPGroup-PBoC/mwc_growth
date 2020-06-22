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
colors, palette = mwc.viz.plotting_style()

force_compile = True

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
model = mwc.bayes.StanModel('../stan/constrained_entropy_estimation.stan', force_compile=force_compile)
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

# %%

# Define constants for the curves
rep_range = np.logspace(0, 3, 200)
temp_range = np.linspace(30, 45, 200) + 273.15
# Set up the figure canvases
fig, ax = plt.subplots(2, 3, figsize=(6, 4))
for a in ax:
    mwc.viz.despine(a)

for i in range(3):
    ax[0, i].set_xscale('log')
    ax[0, i].set_xlabel('repressors per cell')
    ax[0, i].set_yscale('log')
    ax[0, i].set_ylabel('fold-change')
    ax[1, i].set_xlabel('1000 / T [K$^{-1}$]')
    ax[1, i].set_ylabel(r'$\varepsilon$ [$k_BT$]')

temp_colors = {32:colors['blue'], 37: colors['purple'], 42:colors['red']}
for i, t in enumerate(data_dict['temp'] - 273.15):
    _ref = lf_samples[lf_samples['reference_temp']==t]
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
        ax[1, 0].hist(_ref[_ref['parameter']=='delta_SR']['value'].values, 
                      color=temp_colors[g], alpha=0.2, bins=75) 
        for j, r in enumerate(rep_range):
            arch = (1 + (r/4.6E6) * np.exp(-epra))**-1
            cred_region[:,  j] = mwc.stats.compute_hpd(arch, 0.95)
        ax[0, i].fill_between(rep_range, cred_region[0, :], cred_region[1, :],
                            color=temp_colors[g], alpha=0.45, label='__nolegend__')
plt.tight_layout()

    # %%
