# -*- coding: utf-8 -*-
import sys
import os
import numpy as np 
import pandas as pd 
import pystan 
import glob
sys.path.insert(0, '../../../../')
import mwc.bayes
import mwc.stats
import mwc.io
import mwc.process
import mwc.model
import mwc.validation
constants = mwc.model.load_constants()

# Define the experimental constants
DATE = 20190301
RUN_NO = 1
TEMP = 37
CARBON = 'acetate'
OPERATOR = 'O2' 
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
mean_auto = (auto_df['fluor2_mean_death']  ).mean()

# Compute the fluctuations. 
family_df = mwc.process.family_reunion(growth_df, fluo_channel=2)

# Insert identifying information.
family_df['date'] = DATE
family_df['run_no'] = RUN_NO
family_df['temp'] = TEMP
family_df['carbon'] = CARBON
family_df['operator'] = OPERATOR

# Remove cells which do not have any mesured intensity
family_df = family_df[((family_df['I_1'] - mean_auto) * family_df['area_1'] > 0) & ((family_df['I_2'] - mean_auto) * family_df['area_2'] > 0)]

# Save the fluctuations to output. 
family_df.to_csv(f'output/{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{OPERATOR}_fluctuations.csv', index=False) 

# Assemble the data dictionary and sample. 
data_dict = dict(N=len(family_df), 
                I1=(family_df['I_1'] - mean_auto) * family_df['area_1'], 
                I2=(family_df['I_2'] - mean_auto) * family_df['area_2'])

# Perform inference of calibration factor. 
model = mwc.bayes.loadStanModel('../../../stan/calibration_factor.stan')
samples = model.sampling(data_dict, iter=5000, chains=4)
samples_df = samples.to_dataframe()
samples_df = samples_df[['alpha', 'lp__']]
samples_df.rename(columns={'lp__':'log_prob'}, inplace=True)
samples_df.to_csv(f'output/{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{OPERATOR}_cal_factor_samples.csv', index=False)

# Generate the dilution summary figure. 
_ = mwc.validation.dilution_summary(family_df, mean_auto, samples_df['alpha'], ip_dist=IP_DIST, fname='dilution_summary',
                                    title=f'{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{OPERATOR}')

# Compute the fold-change. 
# --------------------------
# Load delta data.
delta_df = mwc.process.parse_clists(glob.glob(f'{data_dir}snaps/delta_00ngml/xy*/clist.mat'))
delta_df = mwc.process.morphological_filter(delta_df, ip_dist=IP_DIST)

# Compute the mean YFP value for autofluorescence and delta LacI
mean_auto_yfp = (auto_df['fluor1_mean_death']).mean()
mean_delta_yfp = (delta_df['fluor1_mean_death']).mean() - mean_auto_yfp

# Compute the mode and hpd of the calibration factor. 
alpha_mode = samples_df.iloc[np.argmax(samples_df['log_prob'].values)]['alpha']
alpha_median = np.median(samples_df['alpha'])
alpha_mean = np.mean(samples_df['alpha'])
hpd_min, hpd_max = mwc.stats.compute_hpd(samples_df['alpha'], 0.95)

# Iterate through all concentrations. 
concs = glob.glob(f'{data_dir}snaps/*')
dfs = []
for i, c in enumerate(concs):
    clists = glob.glob(f'{c}/xy*/clist.mat')
    _df = mwc.process.parse_clists(clists)
    _df = mwc.process.morphological_filter(_df, ip_dist=IP_DIST)
    strain, atc= c.split('/')[-1].split('_')
    atc = float(atc.split('ngml')[0])

    # Extract the important data
    _df['strain'] = strain
    _df['mean_bg_yfp'] = _df['fluor1_bg_death']
    _df['mean_bg_mCherry'] = _df['fluor2_bg_death']
    _df['mean_yfp'] = _df['fluor1_mean_death'] 
    _df['mean_mCherry'] = _df['fluor2_mean_death'] 
    _df['atc_ngml'] = atc
    _df['date'] = DATE
    _df['area_pix'] = _df['area_death']
    _df['carbon'] = CARBON
    _df['fold_change'] = (_df['mean_yfp'] - mean_auto_yfp) / mean_delta_yfp
    _df['temp'] = TEMP
    _df['operator'] = OPERATOR
    _df['run_number'] = RUN_NO
    _df['yfp_bg_val'] = _df['fluor1_bg_death']
    _df['mCherry_bg_val'] = _df['fluor2_bg_death']
    _df['alpha_mode'] = alpha_mode
    _df['alpha_median'] = alpha_median
    _df['alpha_mean'] = alpha_mode
    _df['alpha_hpd_min'] = hpd_min
    _df['alpha_hpd_max'] = hpd_max
    dfs.append(_df[['strain', 'area_pix', 'mean_yfp', 'mean_mCherry', 
                   'atc_ngml', 'date', 'carbon', 'temp', 'operator', 'run_number', 'yfp_bg_val', 'mCherry_bg_val',
                   'alpha_mode', 'alpha_mean', 'alpha_median', 'alpha_hpd_min', 'alpha_hpd_max',
                   'fold_change']])
fc_df = pd.concat(dfs)

# Save to disk. 
fc_df.to_csv(f'output/{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{OPERATOR}_foldchange.csv', index=False)

# Generate the fold-change summary figure
_ = mwc.validation.fc_summary_microscopy(fc_df, samples_df, operator='O2', constants=constants, fname='foldchange_summary', 
                                         title=f'{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{OPERATOR}')
#

