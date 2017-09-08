import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import glob
import sys
sys.path.insert(0, '../../../')
import mwc_growth as mwc
mwc.set_plotting_style()
sns.set_palette('deep', n_colors=8)

# Set the constants.
DATE = '20170814'
RUN_NO = 'r1'


# Load the data
csv_files = glob.glob('output/*titration*.csv')

data = pd.read_csv(csv_files[0], comment='#')

# group the data by IPTG concentration, aTc, and strain.
grouped = data.groupby(['iptg_um'])

dfs = []
for g, d in grouped:
    auto_val = d[d['strain'] == 'auto']['mean_YFP_H'].values
    delta_val = d[d['strain'] == 'delta']['mean_YFP_H'].values
    _d = d.copy()
    _d['foldchange'] = (d['mean_YFP_H'] - auto_val) / (delta_val - auto_val)
    dfs.append(_d)

# Concatenate the data and generate the sanity-check plot.
_data = pd.concat(dfs, axis=0, ignore_index=True)
target_parts = csv_files[0].split('_')[:-3]
target_parts.append('foldchange.csv')
target = "_".join(target_parts)
if os.path.exists(target):
    print('Deleting previous version of file.')
    os.remove(target)
with open('comments.txt', 'r') as f:
    comments = f.readlines()
with open(target, 'a') as f:
    for line in comments:
        f.write(line)
    _data.to_csv(f)


# Generate a plot of auto and delta values as a function of IPTG.
plt.figure()
data_slc = _data[(_data['strain']=='auto') | (_data['strain']=='delta')]
grouped = data_slc.groupby('strain')
for g, d in grouped:
    plt.semilogx(d['iptg_um'], d['mean_YFP_H'], '--o', label=g)
plt.xlabel('IPTG [µM]')
plt.ylabel('mean YFP')
plt.legend(loc='upper left')
plt.savefig('output/{0}_{1}_auto_delta_consistency.png'.format(DATE, RUN_NO), bbox_inches='tight')


data_slc = _data[(_data['strain']!='auto') & (_data['strain']!='delta')]

plt.figure()
grouped = data_slc.groupby(['atc_ngml'])
for g, d in grouped:
    plt.semilogx(d['iptg_um'], d['foldchange'], '-o', lw=1, label=g)
plt.xlabel('IPTG [µM]')
plt.ylabel('fold-change')
plt.legend(loc='upper left', title='aTc (ng/mL)')
plt.ylim([-.1, 1.2])
plt.savefig('output/{0}_{1}_intensity_IPTG_titration.png'.format(DATE, RUN_NO), bbox_inches='tight')
