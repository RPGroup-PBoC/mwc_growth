# -*- coding: utf-8 -*-
import sys
import os
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import pystan 
import glob
sys.path.insert(0, '../../../../')
import mwc.viz
import mwc.bayes
import mwc.stats
import mwc.io
import mwc.process
mwc.viz.personal_style()

# Define the experimental constants
DATE = 20181002
RUN_NO = 1
TEMP = 37
CARBON = 'glucose'
OPERATOR =  'O2'
IP_DIST = 0.065

# Make the output directory if necessary
if os.path.isdir('output') == False:
    os.mkdir('output')

# Define the data directory 
data_dir = f'../../../../data/images/{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{OPERATOR}_dilution/'

# Process the growth files. 
growth_clists = glob.glob(f'{data_dir}growth/xy*/clist.mat')
growth_df = mwc.process.parse_clists(growth_clists)

# Pass through the morphological filter
growth_df = mwc.process.morphological_filter(growth_df, ip_dist= IP_DIST)

# Remove cells with an error frame. 
growth_df['valid'] = growth_df['error_frame'].isnull()
growth_df = growth_df[growth_df['valid'] == True]

# Load the autofluorescence data. 
auto_df = mwc.process.parse_clists(glob.glob(f'{data_dir}snaps/auto_00ngml/xy*/clist.mat'))
auto_df = mwc.process.morphological_filter(auto_df, ip_dist=IP_DIST)
mean_auto = auto_df['fluor2_mean_death'].mean()

# Compute the fluctuations. 
fluct_df = mwc.process.compute_fluctuations(growth_df, auto_val=mean_auto, fluo_key='fluor2_mean_death')

# Insert identifying information.
fluct_df['date'] = DATE
fluct_df['run_no'] = RUN_NO
fluct_df['temp'] = TEMP
fluct_df['carbon'] = CARBON
fluct_df['operator'] = OPERATOR

# Remove cells which do not have any mesured intensity
fluct_df = fluct_df[(fluct_df['I_1'] > 0) & (fluct_df['I_2'] > 0)]

# Save the fluctuations to output. 
fluct_df.to_csv(f'output/{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{OPERATOR}_fluctuations.csv', index=False) 

# Perform inference of calibration factor. 
model = mwc.bayes.loadStanModel('../../../stan/calibration_factor.stan')

# Assemble the data dictionary and sample. 
data_dict = dict(N=len(fluct_df), I1=fluct_df['I_1'], I2=fluct_df['I_2'])
samples = model.sampling(data_dict, iter=5000, chains=4)
samples
samples_df = samples.to_dataframe()
samples_df = samples_df[['alpha', 'lp__']]
samples_df.rename(columns={'lp__':'log_prob'}, inplace=True)
samples_df.to_csv(f'output/{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{OPERATOR}_cal_factor_samples.csv', index=False)

# Compute the fold-change. 
# --------------------------
# Load delta data.
delta_df = mwc.process.parse_clists(glob.glob(f'{data_dir}snaps/delta_00ngml/xy*/clist.mat'))
delta_df = mwc.process.morphological_filter(delta_df, ip_dist=IP_DIST)

# Compute the mean YFP value for autofluorescence and delta LacI
mean_auto_yfp = auto_df['fluor1_mean_death'].mean()
mean_delta_yfp = delta_df['fluor1_mean_death'].mean() - mean_auto_yfp

# Iterate through all concentrations. 
concs = glob.glob(f'{data_dir}snaps/dilution*')
dfs = []
for i, c in enumerate(concs):
    clists = glob.glob(f'{c}/xy*/clist.mat')
    _df = mwc.process.parse_clists(clists)
    _df = mwc.process.morphological_filter(_df, ip_dist=IP_DIST)
    
    # Extract the important data
    _df['mean_yfp'] = _df['fluor1_mean_death'] - mean_auto_yfp
    _df['mean_mCherry'] = _df['fluor2_mean_death'] - mean_auto
    _df['fold_change'] = _df['mean_yfp'] / mean_delta_yfp
    _df['repressors'] = _df['mean_mCherry'] * _df['area_death'] / np.mean(samples_df['alpha'])
    _df['atc_ngml'] = float(c.split('dilution_')[1].split('ngml')[0])
    _df['date'] = DATE
    _df['carbon'] = CARBON
    _df['temp'] = TEMP
    _df['operator'] = OPERATOR
    _df['run_number'] = RUN_NO
    dfs.append(_df[['mean_yfp', 'mean_mCherry', 'fold_change', 'repressors',
                   'atc_ngml', 'date', 'carbon', 'temp', 'operator', 'run_number']])
fc_df = pd.concat(dfs)

# Save to disk. 
fc_df.to_csv(f'output/{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{OPERATOR}_foldchange.csv', index=False)

# Compute the hpd of alpha. 
hpd = mwc.stats.compute_hpd(samples_df['alpha'], 0.95)

# Bin the fluctuation data. 
binned = mwc.stats.bin_by_events(fluct_df, 50)

# Set up the figure canvas. 
fig, ax = plt.subplots(1, 2, figsize=(4, 2))

# Add appropriate labels and formatting. 
ax[0].set_xlabel('$I_1 + I_2$', fontsize=8)
ax[0].set_ylabel('$(I_1 - I_2)^2$', fontsize=8)
ax[1].set_xlabel('calibration factor [a.u. / molecule]', fontsize=8)
ax[1].set_ylabel('probability')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlim([1E0, 5E5])

# Plot the fluctuations.
_ = ax[0].plot(fluct_df['summed'], fluct_df['fluct'], 'k.', ms=0.5, alpha=0.75, label='raw data')
_ = ax[0].plot(binned['summed'], binned['fluct'], '.', color='dodgerblue', ms=3, label='average')

# Plot the credible region
summed_range = np.logspace(0, 6, 200) 
_ = ax[0].fill_between(summed_range, hpd[0] * summed_range, hpd[1] * summed_range, color='tomato', 
                       alpha=0.5, label='credible region', zorder=100)


# Plot the samples. 
_ = ax[1].hist(samples_df['alpha'], bins=100, color='slategray', edgecolor='k', lw=0.05)
_ = ax[1].fill_betweenx(np.linspace(0, ax[1].get_ylim()[1]), hpd[0], hpd[1], color='tomato', alpha=0.25, zorder=100)

# Add legend, format axes, and save. 
ax[0].legend(fontsize=6)
mwc.viz.format_axes()
plt.tight_layout() 
plt.savefig(f'output/{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{OPERATOR}_dilution.png', bbox_inches='tight')

