# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import scipy.stats
import sys
sys.path.insert(0, '../../../')
import mwc.viz
import mwc.process
import mwc.stats
import mwc.bayes
import seaborn as sns
import imp
mwc.viz.personal_style()
%matplotlib inline
imp.reload(mwc.viz)
# Define the experimental parameters.
DATE = 20180209
TEMP = 37  # in °C
CARBON = 'glucose'
OPERATOR = 'O2'
MICROSCOPE = 'tenjin'

# ############################
# Nothing below here should change
# ############################

# %% Processing of data
data_dir = '../../../data/images/{}_{}_{}C_{}_{}_dilution/'.format(
    DATE, MICROSCOPE, TEMP, CARBON, OPERATOR)

# Extract file names and parse.
growth_files = glob.glob('{}growth*/xy*/clist.mat'.format(data_dir))
excluded_props = ['Fluor2 mean death']
growth_df = mwc.process.parse_clists(
    growth_files, excluded_props=excluded_props)


# %%
imp.reload(mwc.process)
snap_groups = glob.glob('{}/snaps*'.format(data_dir))
excluded_props = ['Area death', 'Cell ID', 'Cell birth time', 'Cell death time',
                  'Daughter1 ID', 'Daughter2 ID', 'Mother ID']
snap_groups
snap_dfs = []
for i, s in enumerate(snap_groups):
    _, strain, atc_conc = s.split('/')[-1].split('_')
    atc_conc = float(atc_conc.split('ngmL')[0])
    added_props = {'strain': strain, 'atc_conc_ngmL': atc_conc}
    clists = glob.glob('{}/xy*/clist.mat'.format(s))
    _df = mwc.process.parse_clists(clists, added_props=added_props,
                                   excluded_props=excluded_props)
    snap_dfs.append(_df)
snap_df = pd.concat(snap_dfs, ignore_index=True)

# %% Computation of fluctuations.
auto_val = np.mean(snap_df[snap_df['strain'] ==
                           'autofluorescence']['fluor1_mean_death'])
fluct_df = mwc.process.compute_fluctuations(growth_df, auto_val)

# %% Estimate the calibration factor.
alpha_opt, alpha_err = mwc.bayes.estimate_calibration_factor(growth_df['I_1'],
                                                             growth_df['I_2'])


# %% Plot calibration factor summary statistic plot.

# Set the ranges for the plots
min_summed = growth_df['summed'].min()
max_summed = growth_df['summed'].max()
summed_range = np.logspace(min_summed, max_summed, 500)
min_alpha = alpha_opt - 0.75 * alpha_opt
max_alpha = alpha_opt + 0.75 * alpha_opt
alpha_range = np.linspace(min_alpha, max_alpha, 500)

# Compute the gaussian approximation and normalize the log posterior
approx_pdf = scipy.stats.norm.pdf(alpha_range, alpha_opt, alpha_err)

# Evaluate the log posterior over the alpha range.
log_post = np.zeros_like(alpha_range)
for i, a in enumerate(alpha_range):
    log_post[i] = mwc.bayes.deterministic_log_posterior(
        a, growth_df['I_1'], growth_df['I_2'], neg=False)
post = np.exp(log_post - np.sum(log_post))

# Bin the growth data for a sanity check.
bin_size = 50  # arbitrary choice.
avg_sum, avg_fluc = mwc.stats.bin_by_events(growth_df, bin_size)

# Set up the figure canvas.
fig, ax = plt.subplots(1, 2, figsize=(6, 3))
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlabel(r'$I_1 + I_2$ [a.u.]')
ax[0].set_ylabel(r'$(I_1 - I_2)^2$ [a.u.]')
ax[1].set_xlabel(r'$\alpha$ [a.u.] / mol.')
ax[1].set_ylabel(r'$\propto f(\alpha)$')

# Plot the fluctuations
_ = ax[0].plot(growth_df['summed'], growth_df['fluct'], '.', color='slategray',
               alpha=0.4, label='data')
_ = ax[0].plot(avg_sum, avg_fluct, 'o', label='binned data')
_ = ax[0].plot(summed_range, alpha_opt * summed_range, color='dodgerblue',
               label='$\alpha$ = {0:0f} $\pm$ {1:0f} a.u. / mol.'.format(alpha_opt, alpha_err))
_ = ax[0].legend(loc='upper left')

# Plot the posterior and approximation.
_ = ax[1].plot(alpha_range, post, color='dodgerblue', label='posterior')
_ = ax[1].fill_between(alpha_range, post, color='dodgerblue', alpha=0.4,
                       label='__nolegend__')
_ = ax[1].plot(alpha_range, approx_pdf, '--', lw=1,
               color='tomato', label='Gaussian approx.')
_ = ax[1].legend(loc='upper right')

# Format and save
sns.despine(offset=5, trim=True)
plt.tight_layout()
plt.savefig('/Users/gchure/Desktop/test.pdf', bbox_inches='tight',
            transparent=True)
