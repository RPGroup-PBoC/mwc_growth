# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pystan
import glob
sys.path.insert(0, '../../../../')
import mwc.bayes
import mwc.viz
import mwc.stats
import matplotlib.gridspec as gridspec
mwc.viz.personal_style()

# Define the experimental constants.
DATE = 20180908
RUN = 'r1'
TEMP = 37  # in °C
CARBON = 'glucose'
OPERATOR = 'O2'

BASE_NAME = '{}_{}_{}C_{}_{}'.format(DATE, RUN, TEMP, CARBON, OPERATOR)

# %%
# ################################
# Nothing below here should change
# ################################
if os.path.isdir('output') == False:
    os.mkdir('output')

# Load and trim the data 
data = pd.read_csv('{}_growth.csv'.format(BASE_NAME))
dfs = []
for g, d in data.groupby('sample'):
    d.sort_values('growth_time', inplace=True)
    d['rel_abs'] = d['od_600nm'].values / d.iloc[0]['od_600nm']
    dfs.append(d)
data = pd.concat(dfs)

#%% Load the stan model. 
model = pystan.StanModel('../../../stan/growth_rate.stan')

# %% Sample
data_dict = dict(J=len(data['sample'].unique()), N=len(data), idx=data['sample'].values.astype(int),
            time=data['growth_time'], rel_abs=data['rel_abs'])
samples = model.sampling(data_dict, iter=8000, chains=4, pars=['lambda', 'sigma', 
            'lambda_mu', 'lambda_sigma'])
print('finished!')
samples
# %%
with pm.Model() as model:
    a0 = pm.HalfNormal('a0', sd=1)
    lam = pm.HalfFlat('lambda')
    gamma = mwc.bayes.Jeffreys('gamma', lower=1E-9, upper=100)

    # Compute the expected value.
    time = data['elapsed_time_min'].values
    mu = np.log(a0) + time * lam

    # Define the likelihood and sample.
    like = pm.Cauchy('like', mu, gamma, observed=np.log(
        data['OD_600'].values))
    trace = pm.sample(tune=5000, draws=5000)
    trace_df = mwc.stats.trace_to_dataframe(trace, model)
    stats = mwc.stats.compute_statistics(trace_df)

# %% Compute the best fit and credible region
modes = {}
hpds = {}
grouped = stats.groupby('parameter')
for g, d in grouped:
    modes[g] = d['mode'].values[0]
    hpds[g] = [d[['hpd_min', 'hpd_max']].values]

time_range = np.linspace(data['elapsed_time_min'].min(),
                         data['elapsed_time_min'].max(), 500)
best_fit = modes['a0'] * np.exp(modes['lambda'] * time_range)

cred_region = np.zeros((2, len(time_range)))
for i, t in enumerate(time_range):
    theo = trace_df['a0'] * np.exp(t * trace_df['lambda'])
    cred_region[:, i] = mwc.stats.compute_hpd(theo, mass_frac=0.95)

# %% Compute the doubling time.
t_double = np.log(2) / modes['lambda']
t_double_hpd = mwc.stats.compute_hpd(np.log(2) / trace_df['lambda'],
                                     mass_frac=0.95)

#%% Plot the fit.
fig = plt.figure(figsize=(5, 3))
gs = gridspec.GridSpec(2, 4)
ax0 = plt.subplot(gs[:2, :2])
ax1 = plt.subplot(gs[0, 2:])
ax2 = plt.subplot(gs[1, 2:])
fig.text(0.13, 0.92, '$t_\mathrm{double} = %d_{-%d}^{ + %d}$ min'
         % (t_double, t_double - t_double_hpd[0], t_double_hpd[1] - t_double),
         fontsize=8)
ax1.set_title('sampling frequency')
ax0.set_xlabel('time [min]')
ax0.set_ylabel('OD_600')
ax1.set_xlabel(r'$\lambda^{-1}$ [min]')
ax2.set_xlabel('$A_0$ [a.u.]')
_ = ax0.plot(data['elapsed_time_min'], data['OD_600'], '.',
             ms=2, label='data')
_ = ax0.plot(time_range, best_fit, label='fit', lw=1, alpha=0.8)
_ = ax0.fill_between(time_range, cred_region[0, :], cred_region[1, :],
                     alpha=0.3, color='firebrick', label='__nolegend__')
_ = ax0.legend()

# Plot the sampling distributions.
_ = ax1.hist(1 / trace_df['lambda'], bins=80, normed=True, alpha=0.5,
             lw=1, edgecolor='k', histtype='stepfilled')
_ = ax2.hist(trace_df['a0'], bins=80, normed=True, alpha=0.5, lw=1,
             edgecolor='k', histtype='stepfilled')

plt.tight_layout()
sns.despine(offset=7)
plt.savefig('output/{}_growth.png'.format(BASE_NAME),
            bbox_inches='tight')
