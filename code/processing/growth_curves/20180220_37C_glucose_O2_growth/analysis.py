# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pymc3 as pm
import glob
import sys
import os
sys.path.insert(0, '../../../../')
import mwc.bayes
import mwc.viz
import mwc.stats
import matplotlib.gridspec as gridspec
mwc.viz.personal_style()

# Define the experimental constants.
DATE = 20180220
TEMP = 37  # in °C
CARBON = 'glucose'
OPERATOR = 'O2'

# ################################
# Nothing below here should change
# ################################
if os.path.isdir('output') == False:
    os.mkdir('output')


# Load and trim the data to start at 0D = 0.05
data = pd.read_csv('{}_{}C_{}_{}_growth.csv'.format(DATE, TEMP, CARBON,
                                                    OPERATOR))
data = data[data['absorbance'] >= 0.05]

with pm.Model() as model:
    a0 = pm.HalfNormal('a0', sd=1)
    lam = pm.HalfFlat('lambda')
    gamma = mwc.bayes.Jeffreys('gamma', lower=1E-9, upper=100)

    # Compute the expected value.
    time = data['elapsed_time_min'].values
    mu = np.log(a0) + time * lam

    # Define the likelihood and sample.
    like = pm.Cauchy('like', mu, gamma, observed=np.log(
        data['absorbance'].values))
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
ax0.set_ylabel('absorbance')
ax1.set_xlabel(r'$\lambda^{-1}$ [min]')
ax2.set_xlabel('$A_0$ [a.u.]')
_ = ax0.plot(data['elapsed_time_min'], data['absorbance'], '.',
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
plt.savefig('output/{}_{}C_{}_{}_growth.png'.format(DATE,
                                                    TEMP, CARBON, OPERATOR), bbox_inches='tight')
