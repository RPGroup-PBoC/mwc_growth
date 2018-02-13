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
import mwc.model
import seaborn as sns
%matplotlib inline

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
mcherry_auto_val = np.mean(snap_df[snap_df['strain'] ==
                                   'autofluorescence']['fluor1_mean_death'])

yfp_auto_val = np.mean(snap_df[snap_df['strain'] ==
                               'autofluorescence']['fluor2_mean_death'])
fluct_df = mwc.process.compute_fluctuations(growth_df, mcherry_auto_val)
fluct_df.to_csv('output/{}_{}_{}C_{}_{}_fluctuations.csv'.format(DATE,
                                                                 MICROSCOPE, TEMP, CARBON, OPERATOR))

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
plt.savefig('output/{}_{}_{}C_{}_{}_calibration_factor.png'.format(DATE, MICROSCOPE, TEMP,
                                                                   CARBON, OPERATOR),
            bbox_inches='tight', transparent=True)

# %% Compute the fold-change for the other samples.
# Subtract the autofluorescence from the snap dataframe.
snap_df['fluor1_sub'] = snap_df['death_area'] * \
    (snap_df['fluor1_mean_death'] - mcherry_auto_val)
snap_df['fluor2_sub'] = snap_df['death_area'] * \
    (snap_df['fluor2_mean_death'] - yfp_auto_val)

# Compute the mean expression for ΔLacI.
mean_delta_yfp = snap_df[snap_df['strain'] == 'deltaLacI']['fluor2_sub'].mean()

# Group the dilution strains by their atc concentration.
dilution = snap_df[snap_df['strain'] == 'dilution']
grouped = dilution.groupby('atc_conc_ngmL')

# Set up the fold-change dataframe.
fc_df = pd.DataFrame([], columns=['atc_conc_ngmL',
                                  'mean_repressors', 'mean_yfp', 'fold_change'])

# Loop through each concentration and compute the fold-change.
for g, d in grouped:
    mean_rep = np.mean(d['fluor1_sub'] / alpha_opt)
    mean_yfp = np.mean(d['fluor2_sub'])
    fc = mean_yfp / mean_delta_yfp
    conc_dict = {'atc_conc_ngmL': g, 'mean_repressors': mean_rep, 'mean_yfp': mean_yfp,
                 'fold_change': fc}
    fc_df = fc.append(conc_dict, ignore_index=True)

# Save the fold-change Dataframe to output.
fc_df.to_csv('output/{}_{}_{}C_{}_{}_foldchange.csv'.format(DATE,
                                                            MICROSCOPE, TEMP, CARBON, OPERATOR),
             index=False)

# %% Plot the fold-change curve for a sanity check.
OP_EN = {'O1': -15.0, 'O2': -13.9, 'O3': -9.3}  # in units of kBT.

# Set up the architectural parameters.
c = 0  # Allosteric effector concentration.
ka = 139
ki = 0.53
ep_ai = 4.5  # in kBT
rep_range = np.logspace(0, 4, 500)
ep_r = OP_EN[OPERATOR]  # in kBT.

# Compute the theoretical fold-change.
arch = mwc.model.SimpleRepression(
    rep_range, ep_r, effector_conc=c, ka=ka, ki=ki, ep_ai=ep_ai)
theo_fc = arch.fold_change()

# Set up the figure and plot the result.
fig, ax = plt.subplots(1, 1, figsize=(4, 3))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number of repressors')
ax.set_ylabel('fold-change')
ax.set_ylim([1E-4, 1.2])
# Plot the theory and data.
_ = ax.plot(rep_range, theo_fc, label='prediction')
_ = ax.plot(fc_df['mean_repressors'], fc_df['fold_change'], 'o', label='data')
_ = ax.legend()
mwc.viz.format_axes()
plt.tight_layout()
plt.savefig('output/{}_{}_{}C_{}_{}_foldchange.png'.format(DATE, MICROSCOPE,
                                                           TEMP, CARBON, OPERATOR), bbox_inches='tight', transparent=True)
