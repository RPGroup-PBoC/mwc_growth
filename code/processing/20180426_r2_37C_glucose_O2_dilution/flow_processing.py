# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import fcsparser
import pandas as pd
import os
import glob
import sys
sys.path.insert(0, '../../../')
import mwc.viz
import mwc.flow
import mwc.validation
import scipy.stats
colors = mwc.viz.personal_style()

# Define the details fo the expriment.
DATE = 20180426
CARBON = 'glucose'
TEMP = 37  # in C
OPERATOR = 'O2'


# -----------------
# %% Renaming
# ------------------
FCS_PATTERN = 'RP{}-{}-{}'.format(str(DATE)[:4], str(DATE)[4:6], str(DATE)[6:])
ATC_CONC = [0, 0, 1, 2, 3, 4, 7, 10]  # in ng/mL
IPTG_CONC = [0, 0.1, 5, 10, 25, 50, 75, 100, 250, 500, 1000, 5000]  # in µM

savedir = '../../../data/flow/'

# Define the order of rows and the cols.
ROWS = ('auto', 'delta', 'dilution', 'dilution',
        'dilution', 'dilution', 'dilution', 'dilution')
# Get the names of the files.
files = glob.glob('../../../data/flow/{}*.fcs'.format(FCS_PATTERN))
files = np.sort(files)

# Break the list up into columns.
ncols, nrows = len(IPTG_CONC), len(ROWS)
col_groups = [files[i:i + nrows] for i in range(0, len(files), nrows)]
for i, col in enumerate(col_groups):
    for j, samp in enumerate(col):
        # Define the new name.
        name = '{}_{}C_{}_{}_{}_{}ngmL_atc_{}uM_iptg'.format(DATE, TEMP, CARBON, OPERATOR, ROWS[j],
                                                             ATC_CONC[j], IPTG_CONC[i])

        # Load the file using fcsparser and save to csv.
        _ = mwc.flow.fcs_to_csv(samp, '{}{}.csv'.format(savedir, name))

# ------------------
# %% Processing
# ------------------
# Validate the flow files.
csv_files = glob.glob('../../../data/flow/{}*.csv'.format(DATE))
_ = mwc.validation.FlowValidation(csv_files).run_suite()

# Load all data and assemble into one dataframe.
df = pd.DataFrame([], columns=['date', 'operator', 'carbon', 'temp', 'strain',
                               'atc_ngml', 'iptg_um', 'mean_fitc'])
for i, f in enumerate(csv_files):
    # Get the information.
    _, _, _, _, strain, atc, _,  iptg, _ = f.split('/')[-1].split('_')
    atc = float(atc[:-4])
    iptg = float(iptg[:-2])

    # Load the dataframe and prune.
    _df = pd.read_csv(f)
    gated = _df[_df['gate'] == 1]['FITC-A'].mean()
    df = df.append({'date': DATE, 'operator': OPERATOR, 'carbon': CARBON, 'temp': TEMP,
                    'strain': strain, 'atc_ngml': atc, 'iptg_um': iptg, 'mean_fitc': gated},
                   ignore_index=True)

# Compute the fold-change
mean_auto = df[df['strain'] == 'auto'].groupby(
    ['iptg_um'])['mean_fitc'].mean()
denom = df[df['strain'] == 'delta'].groupby(
    ['iptg_um'])['mean_fitc'].mean() - mean_auto
denom = denom.to_dict()
mean_auto = mean_auto.to_dict()
mean_auto = [mean_auto[z] for z in df['iptg_um'].values]
denom = [denom[z] for z in df['iptg_um'].values]
df['fold_change'] = (df['mean_fitc'] - mean_auto) / denom

if os.path.exists('./output') == False:
    os.mkdir('./output')
df.to_csv('./output/{}_{}C_{}_{}_flow_events.csv'.format(DATE, TEMP, CARBON, OPERATOR),
          index=False)

# ---------------
# %% Calculate the fold-change
# ---------------
indices = {'auto': 0, 'delta': 1, 'dilution': None}
ones = np.ones((8, 12))
for i, k in enumerate(indices.keys()):
    if (k == 'auto') | (k == 'delta'):
        ones[indices[k], :] *= df[df['strain'] ==
                                  k].sort_values('iptg_um')['fold_change'].values
    else:
        atc_concs = np.sort(df['atc_ngml'].unique())[1:]
        iptg_concs = np.sort(df['iptg_um'].unique())
        for j, a in enumerate(atc_concs):
            for l, p in enumerate(iptg_concs):
                slc = df[(df['strain'] == 'dilution') & (
                    df['atc_ngml'] == a) & (df['iptg_um'] == p)]
                ones[j + 2, l] *= slc['fold_change']


ex_df = pd.read_csv(csv_files[2])
nogate = ex_df[ex_df['gate'] == 0]
gate = ex_df[ex_df['gate'] == 1]

X, Y = np.mgrid[ex_df['FSC-A'].min():ex_df['FSC-A'].max():100j,
                ex_df['SSC-A'].min():ex_df['SSC-A'].max():100j]
positions = np.vstack([X.ravel(), Y.ravel()])
xy = np.vstack([ex_df['FSC-A'], ex_df['SSC-A']])
kernel = scipy.stats.gaussian_kde((xy))
z = kernel(xy)
idx = z.argsort()
xi = np.logspace(3, 5, 100)
yi = np.logspace(2, 6, 100)
a, b, c = ex_df['FSC-A'].values[idx], ex_df['SSC-A'].values[idx], z[idx]
zi = scipy.interpolate.griddata(
    (a, b), c, (xi[None, :], yi[:, None]), method='cubic')

# %%
fig, ax = plt.subplots(1, 2, figsize=(8, 3))
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlabel('forward scatter [a.u.]')
ax[0].set_ylabel('side scatter [a.u.]')
ax[1].set_xlabel('IPTG [µM]')
ax[1].set_ylabel('ATC [ng/mL]')

_ = ax[0].scatter(a, b, c=c, s=0.5, edgecolor='', marker='.', cmap='viridis')
_ = ax[0].plot(gate['FSC-A'], gate['SSC-A'], ',', color='w', alpha=0.4)
mwc.viz.format_axes()
plate = ax[1].imshow(ones, cmap='viridis', aspect='auto')
_ = plt.colorbar(plate, ax=ax[1], label='fold-change')
ax[1].grid(False)
ax[1].set_xticks(np.arange(0, 12, 1))
ax[1].set_yticks(np.arange(0, 8, 1))
_ = ax[1].set_xticklabels(iptg_concs, fontsize=8, rotation=45)
_ = ax[1].set_yticklabels(['auto | 0', 'ΔlacI | 0 ', 'dilution | 1', '2', '3', '4', '7',
                           '10'], fontsize=8)
plt.tight_layout()
plt.savefig('output/{}_{}C_{}_{}_flow_summary.png'.format(DATE,
            TEMP, CARBON, OPERATOR, bbox_inches='tight'))
