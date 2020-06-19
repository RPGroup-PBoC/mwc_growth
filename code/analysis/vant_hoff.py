#%%
import numpy as np 
import pandas as pd
import mwc.bayes
import mwc.stats
import mwc.process
import scipy.stats
import mwc.viz
import matplotlib.pyplot as plt
colors, palette = mwc.viz.plotting_style()

# Load the data 
data = pd.read_csv('../../data/analyzed_foldchange.csv', comment='#')
data = mwc.process.condition_filter(data, carbon='glucose')
data['repressors'] *= 1.16
data = data.groupby(['temp', 'atc_ngml', 
                     'date', 'run_number'])[
                     ['fold_change', 'repressors']].mean().reset_index()


#%%
# Approach 1 - Complete generative model
vanthoff_model = mwc.bayes.StanModel('../stan/vant_hoff.stan', force_compile=True)

#%%
data.loc[data['temp']==32, 'idx'] = 1
data.loc[data['temp']==37, 'idx'] = 2
data.loc[data['temp']==42, 'idx'] = 3
data_dict = {'J':3, 
            'N':len(data),
            'temp_idx':data['idx'].values.astype(int),
            'temps': data['temp'].unique(),
            'R':data['repressors'],
            'Nns':int(4.6E6),
            'foldchange': data['fold_change'].values}
fit, sample = vanthoff_model.sample(data_dict, iter=5000)
fit

# %%
delH = sample['delH'].values
delS = sample['delS'].values
fig, ax = plt.subplots(1, 1, figsize=(3,3))
ax.errorbar(1000/(37+273), sample['epsilon[2]'].mean(), sample['epsilon[2]'].std(),
        color=colors['purple'], fmt='o', linestyle='none', markeredgecolor=colors['grey'],
        markersize=4, markeredgewidth=0.5)
ax.errorbar(1000/(32+273), sample['epsilon[1]'].mean(), sample['epsilon[1]'].std(),
        color=colors['purple'], fmt='o', linestyle='none', markeredgecolor=colors['grey'],
        markersize=4, markeredgewidth=0.5)
ax.errorbar(1000/(42+273), sample['epsilon[3]'].mean(), sample['epsilon[3]'].std(),
        color=colors['purple'], fmt='o', linestyle='none', markeredgecolor=colors['grey'],
        markersize=4, markeredgewidth=0.5)
temp_range = np.linspace(302, 318, 100)

cred_region = np.zeros((2, len(temp_range)))
for i, _t in enumerate(temp_range):
    pred = -(delH/ _t) + delS 
    cred_region[:, i] = mwc.stats.compute_hpd(pred, 0.95)
plt.fill_between(1000/temp_range, cred_region[0, :], cred_region[1, :], alpha=0.5,
        color=colors['purple'])
ax.set_xlabel(r'$\frac{1000}{T}$ [K$^{-1}$]')  
ax.set_ylabel(r'-$\log\frac{K_D}{1\mu\mathrm{M}}$')  
mwc.viz.titlebox(ax, r'$\Delta H = 10.0 \pm 0.1 \, (\times 10^3) k_B$ ; $\Delta S = 17.9 \pm 0.3\, k_B$',
                color=colors['black'])
plt.savefig('./vant_hoff_inference.pdf')

#%%
fig, ax = plt.subplots(1, 1, figsize=(4,4))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('repressors per cell')
ax.set_ylabel('fold-change')

rep_range = np.logspace(0, 3, 200)
summarized = data.groupby(['atc_ngml', 'temp'])[['repressors', 'fold_change']].agg(('mean', 'sem')).reset_index()
for g, d in summarized.groupby(['temp']):
    if g == 32:
        name = 'epsilon[1]'
        c = colors['blue']
    elif g == 37:
        name = 'epsilon[2]'
        c = colors['purple']
    elif g == 42:
        name = 'epsilon[3]'
        c = colors['red']

    ep_high = sample[name].mean() + sample[name].std()
    ep_low = sample[name].mean() - sample[name].std()
    ep, r = np.meshgrid([ep_high, ep_low], rep_range)
    fc = mwc.model.SimpleRepression(r, ep, ka=139, ki=0.53, ep_ai=10000,
                                    effector_conc=0).fold_change()
    ax.fill_between(rep_range, fc[:, 0], fc[:, 1], color=c, alpha=0.5, label='__nolegend__')
    ax.errorbar(d['repressors']['mean'], d['fold_change']['mean'],
                xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'],
                fmt='o', color=c, linestyle='none', ms=4, markeredgecolor=colors['grey'],
                markeredgewidth=0.5, label=f'{int(g)}° C')

ax.legend(fontsize=6)
plt.savefig('./vant_hoff_foldchange.pdf')
#%%
# Approach II: Piecewise model
# Load the stan model
epsilon_model = mwc.bayes.StanModel('../stan/DNA_binding_energy.stan', force_compile=True)
# linreg_model = mwc.bayes.StanModel('../stan/lin_reg.stan', force_compile=True)

# %%
# Assemble the data dictionary
dfs = []
for g, d in data.groupby(['temp']):
    data_dict = {'N':len(d), 
                 'repressors': d['repressors'].values, 
                 'foldchange': d['fold_change'].values,
                 'Nns':int(4.6E6)}
    _, samples = epsilon_model.sample(data_dict)
    samples['temp'] = g
    dfs.append(samples)
epsilon_samples = pd.concat(dfs, sort=False)


# %%
# Infer the slope
epsilon_means = epsilon_samples.groupby(['temp'])['epRA'].mean().reset_index()
epsilon_hpds = {}
for g, d in epsilon_samples.groupby(['temp']):
    low, high = mwc.stats.compute_hpd(d['epRA'], 0.95)
    epsilon_hpds[g] = {'low':low, 'high':high}


# %%
import matplotlib.pyplot as plt
stats = pd.read_csv('../../data/pooled_entropic_parameter_samples.csv', comment='#')
delta_SR = stats[stats['parameter']=='delta_SR']['value'].unique()
ref_epRA = stats[stats['parameter']=='ref_epRA']['value'].unique()

fig, ax = plt.subplots(1, 1, figsize=(3,3))
t_range = np.linspace(300, 330, 200)
cred_region = np.zeros((2, len(t_range)))
for i, t in enumerate(t_range):
    ep = delta_SR * ((37 + 273.15) - t) + ref_epRA 
    cred_region[:, i] = mwc.stats.compute_hpd(ep, 0.95) 

ax.plot(1000 *(epsilon_means['temp'] + 273.15)**-1, epsilon_means['epRA'], 'o',
        markeredgecolor=colors['grey'], markeredgewidth=0.75, zorder=1000)
for t, v in epsilon_hpds.items():
    ax.vlines(1000 / (t + 273.15), v['low'], v['high'], lw=1, color=colors['purple'])
ax.fill_between(1000/t_range, cred_region[0, :], cred_region[1, :], color=colors['purple'],  alpha=0.5)
ax.set_xlabel('1000 / T [K$^{-1}$]')
ax.set_ylabel('$\epsilon / k_BT$')
ax.set_xlim([3.15, 3.30])
# ax.set_xlim([0])
# %%
