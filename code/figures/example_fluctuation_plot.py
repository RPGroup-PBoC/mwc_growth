# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import mwc.stats
mwc.viz.pub_style()


# Load an example data set. 
fluct_data = pd.read_csv('../../data/compiled_fluctuations.csv')
fluct_data = fluct_data[fluct_data['date']==20181002]

# Compute bins of 50 events each. 
bins = mwc.stats.bin_by_events(fluct_data, average=['summed', 'sq_fluct'], bin_size=30)

# Load an example data set. 
samples = pd.read_csv('../../code/processing/microscopy/20181002_r1_37C_glucose_O2_dilution/output/20181002_r1_37C_glucose_O2_cal_factor_samples.csv')

# Compute the theory line
summed_range = np.logspace(1, 6, 200)
alpha_median = np.median(samples['alpha'])
low, high = mwc.stats.compute_hpd(samples['alpha'], 0.95)


fig, ax = plt.subplots(1, 1, figsize=(6, 4))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([1E1, 1E6])
ax.set_xlabel(r'$I_1 + I_2$', fontsize=12)
ax.set_ylabel(r'$(I_1 - I_2)^2$', fontsize=12)
ax.plot(fluct_data['summed'], fluct_data['sq_fluct'], '.', color='k', ms=2, label='raw data')
ax.plot(bins['summed'], bins['sq_fluct'], 'o', color='tomato', ms=5, label='30 events / bin')
ax.plot(summed_range, alpha_median * summed_range, lw=2, color='firebrick', label='α = 490  ± 30 a.u.'.format(alpha_median, np.std(samples['alpha'])))
ax.fill_between(summed_range, low* summed_range, high * summed_range, color='firebrick', alpha=0.3)
ax.legend(loc='lower right')
plt.savefig('../../figs/example_dilution_cloud.svg', bbox_inches='tight')
