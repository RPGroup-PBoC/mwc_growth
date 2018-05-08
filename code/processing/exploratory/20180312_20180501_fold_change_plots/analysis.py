import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import glob
import sys
sys.path.insert(0, '../../../../')
import mwc
mwc.viz.personal_style()
sns.set_palette('deep', n_colors=8)

# Load the data and assemble into one dataframe.
csv_files = glob.glob('./input/*flow_events.csv')
csv_files = np.sort(csv_files)

df = pd.DataFrame([], columns=['date', 'operator', 'carbon', 'temp', 'strain',
                               'atc_ngml', 'iptg_um', 'mean_fitc', 'fold_change'])

# This currently doesn't account for runs. 
# It works only if there are never two runs with the same carbon source on the same day.

for i, f in enumerate(csv_files):
    _df = pd.read_csv(f)
    df = df.append(_df,ignore_index=True)

df = df.sort_values(['date','carbon','atc_ngml','strain','iptg_um']).reset_index(drop=True)


# Generate a plot of fold-change values as a function of IPTG.

# Create dataframe with only dilution strain entries.
dil = df[(df['strain']!='auto') & (df['strain']!='delta')]

# Find the mean and standard error of fold-change per atc concentration and carbon source, 
# averaged over all experiments (excluding auto and delta).
fc_agg = dil.groupby(['carbon','iptg_um','atc_ngml']).fold_change.agg(['mean', 'sem']).reset_index()
fc_agg.rename(columns = {'mean':'fc_mean','sem':'fc_sem'}, inplace = True)


# Plot the mean fold-change values as a function of IPTG, per atc concentration for glucose.
plt.figure()

grouped = fc_agg[fc_agg['carbon']=='glucose'].groupby(['atc_ngml','carbon'])
for group, data in grouped:
    #plt.semilogx(data['iptg_um'], data['fold_change'], '-o', lw=1, label=group) #no error bars
    plt.errorbar(data['iptg_um'], data['fc_mean'], yerr=data['fc_sem'], marker='o', lw=1, label=group[0])
plt.xscale('log')
plt.xlabel('IPTG [µM]')
plt.ylabel('Fold-Change')
plt.legend(loc='upper left', title='aTc (ng/mL)')
plt.ylim([-.1, 1])
plt.title("Fold Change: Glucose")
plt.savefig('output/20180312_20180501_fold_change_IPTG_titration_Glucose.png', bbox_inches='tight')

# Plot the mean fold-change values as a function of IPTG, per atc concentration for glycerol.
plt.figure()

grouped = fc_agg[fc_agg['carbon']=='glycerol'].groupby(['atc_ngml','carbon'])
for group, data in grouped:
    plt.errorbar(data['iptg_um'], data['fc_mean'], yerr=data['fc_sem'], marker='o', lw=1, label=group[0])
plt.xscale('log')
plt.xlabel('IPTG [µM]')
plt.ylabel('Fold-Change')
plt.legend(loc='upper left', title='aTc (ng/mL)')
plt.ylim([-.1, 1])
plt.title("Fold Change: Glycerol")
plt.savefig('output/20180312_20180501_fold_change_IPTG_titration_Glycerol.png', bbox_inches='tight')


# Find the mean and standard error of YFP intensity (fitc) per atc concentration, carbon source, and strain,
# averaged over all experiments (using only auto and delta).

yfp_agg = df[(df['strain']=='auto')|(df['strain']=='delta')].groupby(
    ['strain','carbon','iptg_um']).mean_fitc.agg(['mean', 'sem']).reset_index()
yfp_agg.rename(columns = {'mean':'yfp_mean','sem':'yfp_sem'}, inplace = True)

# Plot mean auto and delta YFP values as a function of IPTG, per strain and carbon source.
plt.figure()

grouped = yfp_agg[yfp_agg['carbon']=='glucose'].groupby('strain')
for group, data in grouped:
    plt.errorbar(data['iptg_um'], data['yfp_mean'], yerr=data['yfp_sem'], marker='o', lw=1, label=group)
plt.xscale('log')
plt.xlabel('IPTG [µM]')
plt.ylabel('FITC')
plt.legend(loc='upper left', title='strain')
plt.ylim([-.1, 50000])
plt.title("YFP: Glucose")
plt.savefig('output/20180312_20180501_YFP_Controls_IPTG_titration_Glucose.png', bbox_inches='tight')

# Plot mean auto and delta YFP values as a function of IPTG, per atc concentration for glycerol.
plt.figure()

grouped = yfp_agg[yfp_agg['carbon']=='glycerol'].groupby('strain')
for group, data in grouped:
    plt.errorbar(data['iptg_um'], data['yfp_mean'], yerr=data['yfp_sem'], marker='o', lw=1, label=group)
plt.xscale('log')
plt.xlabel('IPTG [µM]')
plt.ylabel('FITC')
plt.legend(loc='upper left', title='strain')
plt.ylim([-.1, 50000])
plt.title("YFP: Glycerol")
plt.savefig('output/20180312_20180501_YFP_Controls_IPTG_titration_Glycerol.png', bbox_inches='tight')