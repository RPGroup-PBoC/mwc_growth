#%%
import numpy as np 
import pandas as pd
import mwc.bayes
import mwc.stats
import mwc.process
import scipy.stats

# Load the data 
data = pd.read_csv('../../data/analyzed_foldchange.csv', comment='#')
data = mwc.process.condition_filter(data, carbon='glucose')
data['repressors'] *= 1.16
data = data.groupby(['temp', 'atc_ngml', 
                     'date', 'run_number'])[
                     ['fold_change', 'repressors']].mean().reset_index()


#%%
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
fit, sample = vanthoff_model.sample(data_dict)
fit

# %%
# Load the stan model
epsilon_model = mwc.bayes.StanModel('../stan/DNA_binding_energy.stan')
linreg_model = mwc.bayes.StanModel('../stan/lin_reg.stan', force_compile=True)

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
# data_dict = {'N':3,
#              'mean_vals':epsilon_means['epRA'].values,
#              'temps':epsilon_means['temp'].values + 273.15}
# fit, samples = linreg_model.sample(data_dict)
_out = scipy.stats.linregress(1 / (epsilon_means['temp'] + 273.15),
                                                epsilon_means['epRA'])
delH = _out[0]
delS = _out[1]

# %%
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1, figsize=(3,3))
t_range = np.linspace(300, 330, 200)
ax.plot((epsilon_means['temp'] + 273.15)**-1, epsilon_means['epRA'], 'o')
ax.plot(1/t_range, (delH / t_range) + delS, 'k-')

#%%
# ax.errorbar(t, np.array([ep1['mean'], ep2['mean'], ep3['mean']]), 
            #    np.array([ep1['std'], ep2['std'], ep3['std']]), color=colors['black'],
            #    fmt='o', markeredgecolor=colors['grey'], markeredgewidth=0.5, lw=0.74)
# slope = samples['slope']
# intercept = samples['intercept']
cred_region = np.zeros((2, len(t_range)))
for i, _t in enumerate(t_range):
    pred = -(samples['slope'] / _t) + samples['intercept']
    cred_region[:, i] = mwc.stats.compute_hpd(pred, 0.95)
plt.fill_between(1/t_range, cred_region[0, :], cred_region[1, :], alpha=0.5)
  
# %%
samples


# %%
