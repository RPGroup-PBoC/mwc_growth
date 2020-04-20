#%%
import numpy as np 
import pandas as pd
import mwc.bayes
import mwc.stats
import mwc.process

# Load the data 
data = pd.read_csv('../../data/analyzed_foldchange.csv', comment='#')
data = mwc.process.condition_filter(data, carbon='glucose')
data['repressors'] *= 1.16
data = data.groupby(['temp', 'atc_ngml', 
                     'date', 'run_number'])[
                     ['fold_change', 'repressors']].mean().reset_index()

# %%
# Load the stan model
model = mwc.bayes.StanModel('../stan/vant_hoff.stan', force_compile=True)

# %%
# Assemble the data dictionary
data.loc[data['temp']==32, 'idx'] = 1
data.loc[data['temp']==37, 'idx'] = 2
data.loc[data['temp']==42, 'idx'] = 3
data_dict = {'J': 3,
             'N':len(data), 
             'temps': data['temp'].unique(),
             'temp_idx': data['idx'].values.astype(int),
             'R': data['repressors'].values, 
             'foldchange': data['fold_change'].values,
             'Nns':int(4.6E6)}
fit, samples = model.sample(data_dict)

# %%
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1, figsize=(3,3))
t_range = np.linspace(300, 330, 200)
t = np.array([(273.15 + 32)**-1, (273.15 + 37)**-1, (273.15 + 42)**-1])
ax.plot(t, eps, 'ko')
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
