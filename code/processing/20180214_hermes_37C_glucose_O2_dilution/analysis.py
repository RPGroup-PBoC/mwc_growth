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
import pymc3 as pm
import os
import seaborn as sns
mwc.viz.personal_style()
# Define the experimental parameters.
DATE = 20180214
TEMP = 37  # in °C
CARBON = 'glucose'
OPERATOR = 'O2'
MICROSCOPE = 'hermes'

# ############################
# Nothing below here should change
# ############################
IP_DIST = 0.063
if os.path.exists('./output') == False:
    os.mkdir('./output')

# %% Processing of data
data_dir = '../../../data/images/{}_{}_{}C_{}_{}_dilution/'.format(
    DATE, MICROSCOPE, TEMP, CARBON, OPERATOR)

# Extract file names and parse.
growth_files = glob.glob('{}growth*/xy*/clist.mat'.format(data_dir))
excluded_props = ['Fluor2 mean death']
growth_df = mwc.process.parse_clists(
    growth_files, excluded_props=excluded_props)

# Apply a filter.
growth_df = mwc.process.morphological_filter(growth_df, IP_DIST)

# %% Look at the aspect ratio.
snap_groups = glob.glob('{}/snaps*'.format(data_dir))
excluded_props = ['Area birth', 'Cell ID', 'Cell birth time', 'Cell death time',
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

    print(added_props, len(_df))

    snap_dfs.append(_df)
snap_df = pd.concat(snap_dfs, ignore_index=True)

# Apply area bounds.
snap_df = mwc.process.morphological_filter(snap_df, IP_DIST)

# %% Computation of fluctuations.
auto_strain = snap_df[snap_df['strain'] == 'autofluorescence']
mcherry_auto_val = np.mean(
    auto_strain['fluor1_mean_death'])
yfp_auto_val = np.mean(snap_df[snap_df['strain'] ==
                               'autofluorescence']['fluor2_mean_death'])

fluct_df = mwc.process.compute_fluctuations(growth_df, mcherry_auto_val)
fluct_df.to_csv('output/{}_{}_{}C_{}_{}_fluctuations.csv'.format(DATE,
                                                                 MICROSCOPE, TEMP,
                                                                 CARBON, OPERATOR))
# %% Estimate the calibration factor.
with pm.Model() as model:
    like = mwc.bayes.DeterminsticCalibrationFactor('alpha', I_1=fluct_df['I_1'].values,
                                                   I_2=fluct_df['I_2'].values, testval=1)
    trace = pm.sample(draws=5000, tune=5000, njobs=4)
    trace_df = mwc.stats.trace_to_dataframe(trace, model)
    stats = mwc.stats.compute_statistics(trace_df)
alpha_opt = stats['mode'].values
alpha_max = stats['hpd_max'].values
alpha_min = stats['hpd_min'].values

# %% Plot calibration factor summary statistic plot.
# Set the ranges for the plots
min_summed = np.log10(fluct_df['summed'].min())
max_summed = np.log10(fluct_df['summed'].max())
summed_range = np.logspace(min_summed, max_summed, 500)
min_alpha = alpha_opt - alpha_opt * 0.5
max_alpha = alpha_opt + alpha_opt * 0.5

# Bin the growth data for a sanity check.
bin_size = 50  # arbitrary choice.
avg = mwc.stats.bin_by_events(fluct_df, bin_size)

# Set up the figure canvas.
fig, ax = plt.subplots(1, 2, figsize=(6, 3))
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlabel(r'$I_1 + I_2$ [a.u.]')
ax[0].set_ylabel(r'$(I_1 - I_2)^2$ [a.u.]')
ax[1].set_xlabel(r'$\alpha$ [a.u.] / mol.')
ax[1].set_ylabel('frequency')
ax[0].set_title(r'$\alpha = %d^{+%d}_{-%d}$ a.u. /mol' %
                (alpha_opt[0], alpha_max[0] - alpha_opt[0], alpha_opt[0] - alpha_min[0]))
ax[1].set_title('sampling frequency')


# Plot the fluctuations
_ = ax[0].plot(fluct_df['summed'], fluct_df['fluct'],
               '.', ms=1, alpha=0.4, label='data')
_ = ax[0].plot(avg['summed'], avg['fluct'], '.', label='binned data')
_ = ax[0].plot(summed_range, alpha_opt * summed_range, label='fit')
_ = ax[0].legend()

# Plot the posterior and approximation.
_ = ax[1].hist(trace_df['alpha'], bins=75, alpha=0.5,
               normed=True)
# Plot the mode and HPD ontop of the histogram
ylim = ax[1].get_ylim()[1] / 2
_ = ax[1].plot([alpha_min, alpha_max], [ylim, ylim], color='tomato')
_ = ax[1].plot(alpha_opt, ylim, 'o', color='tomato')

# Format and save
ax[1].set_xlim([min_alpha, max_alpha])
ax[0].set_xlim([100, 10**max_summed])
plt.tight_layout()
mwc.viz.format_axes()
plt.savefig('output/{}_{}_{}C_{}_{}_calibration_factor.png'.format(DATE, MICROSCOPE, TEMP,
                                                                   CARBON, OPERATOR),
            bbox_inches='tight')

# %% Compute the fold-change for the other samples.
# Subtract the autofluorescence from the snap dataframe.
snap_df['fluor1_sub'] = snap_df['area_death'] * \
    (snap_df['fluor1_mean_death'] - mcherry_auto_val)

snap_df['fluor2_sub'] = snap_df['area_death'] * \
    (snap_df['fluor2_mean_death'] - yfp_auto_val)

# Compute the mean expression for ΔLacI.
mean_delta_yfp = snap_df[snap_df['strain'] == 'deltaLacI']['fluor2_sub'].mean()

# Group the dilution strains by their atc concentration.
dilution = snap_df[snap_df['strain'] == 'dilution']
grouped = dilution.groupby('atc_conc_ngml')

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
    fc_df = fc_df.append(conc_dict, ignore_index=True)

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
_ = ax.plot(fc_df['mean_repressors'],
            fc_df['fold_change'], '.', label='data')
_ = ax.legend()
mwc.viz.format_axes()
plt.tight_layout()
plt.savefig('output/{}_{}_{}C_{}_{}_foldchange.png'.format(DATE, MICROSCOPE,
                                                           TEMP, CARBON, OPERATOR),
            bbox_inches='tight')
