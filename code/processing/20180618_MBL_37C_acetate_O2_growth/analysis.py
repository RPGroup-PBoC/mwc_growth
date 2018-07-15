# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import pystan
import matplotlib.pyplot as plt
import glob
sys.path.insert(0,'../../../')
import mwc.viz
import mwc.bayes
colors = mwc.viz.personal_style()

DATE = '20180618'
RUN_NO = 'MBL'
TEMPERATURE = 37 # in C
CARBON = 'glucose'
OPERATOR = 'O2'

# load the data.
data = pd.read_csv(glob.glob('output/*.csv')[0])

#%% Compile the pystan model 
model = pystan.StanModel('../../stan/hierarchical_growth_curve_single_condition.stan')
#%% Set up the data dictionary and sample.
data_dict = {'J':data['colony_idx'].astype(int).max(), 'N':len(data), 
            'colony':data['colony_idx'].astype(int), 'fractional_area':data['fractional_area'],
            'time':data['time_min']}
print('beginning sampling....')
samples = model.sampling(data_dict, iter=5000, chains=4)
print('finished!')
samples
# %% 
sample_df = mwc.bayes.chains_to_dataframe(samples)
stats = mwc.stats.compute_statistics(sample_df)
stats

# Compute the doubling time.
t_double_mode = np.log(2) / stats[stats['parameter']=='lambda_mu']['mode'].values[0]
t_double_max = np.log(2) / stats[stats['parameter']=='lambda_mu']['hpd_min'].values[0]
t_double_min = np.log(2) / stats[stats['parameter']=='lambda_mu']['hpd_max'].values[0]

# Compute the best fit and credible regions from the mode/HPD
time_range = np.linspace(0, data['time_min'].max() + 10, 500)
best_fit = np.exp(time_range * stats[stats['parameter']=='lambda_mu']['mode'].values[0])
upper = np.exp(time_range * stats[stats['parameter']=='lambda_mu']['hpd_max'].values[0])
lower = np.exp(time_range * stats[stats['parameter']=='lambda_mu']['hpd_min'].values[0])

# Set up the plot.
fig, ax = plt.subplots(1, 2, figsize=(6, 3))
ax[0].set_xlabel('time [min]')
ax[0].set_ylabel('fractional area')
ax[1].set_xlabel('$\lambda$ [min$^{-1}$]')
ax[1].set_ylabel(r'$\propto$ probability')


# Plot the best-fit
_ = ax[0].plot(time_range, best_fit, lw=1, color=colors[1], label='fit')
_ = ax[0].fill_between(time_range, lower, upper,  color=colors[1], alpha=0.5, label='__nolegend__')

# Plot the data.

_ = ax[0].plot(data['time_min'], data['fractional_area'], '.', ms=1, color=colors[0], label='microcolonies')
grouped = pd.DataFrame(data.groupby('time_min').apply(mwc.stats.compute_mean_sem, key='fractional_area')).reset_index()
_ = ax[0].errorbar(grouped['time_min'], grouped['mean'], grouped['sem'], lw=1, fmt='.', color=colors[2], label='mean $\pm$ sem', 
markersize=4, zorder=1000)

# Plot the posterior distribution.
for i in range(data['colony_idx'].astype(int).max()):
    if i == 0:
        label = 'λ'
    else:
        label = '__nolegend__'
    samps = sample_df['lambda.{}'.format(i+1)]
    _ = ax[1].plot(np.sort(samps), np.arange(0, len(samps)) / len(samps), color=colors[0], lw=1,
    alpha=0.5, label=label)

_ = ax[1].plot(np.sort(sample_df['lambda_mu']), np.arange(0, len(samps)) / len(samps), color=colors[1],
label='λ$^*$')
ax[1].set_xlim([0, 0.04])
ax[0].legend()

leg = ax[1].legend(title=r'$t_\mathrm{double} = %0.0f^{+%.0f}_{-%.0f}$ min' %(t_double_mode, 
                    t_double_max - t_double_mode, 
                    t_double_mode - t_double_min))
leg.get_title().set_fontsize(8)
plt.tight_layout()
plt.savefig('output/{}_{}_{}C_{}_{}_growth_curves.png'.format(DATE, RUN_NO, TEMPERATURE, CARBON,
OPERATOR))






