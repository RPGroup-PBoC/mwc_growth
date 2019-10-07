#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import mwc.stats
colors, color_list = mwc.viz.personal_style()

# Define the repressor copy number range. 
rep_range = np.logspace(0, 3.75, 200)

# Define the target points for the explanation
example_reps = np.array([10, 30, 100, 500, 1000])
example_facecolors = {10:colors['light_purple'], 30: colors['light_orange'],
                    100:'white', 500:colors['light_red'], 1000:colors['light_blue']}
example_edgecolors = {10:colors['dark_purple'], 30: colors['dark_orange'],
                    100:colors['black'], 500:colors['dark_red'], 1000:colors['dark_blue']}
ref_rep = 100

# Define thermodynamic constants
ep_ai = 4.5 # in kT
ep_RA = -13.9 # in kT
Nns = 4.6E6 # in bp

# Define the fold-change and the deltaF theory curves
pact = (1 + np.exp(-ep_ai))**-1
fc = (1 + pact * (rep_range / Nns) * np.exp(-ep_RA))**-1 
example_fc = (1 + pact * (example_reps/Nns) * np.exp(-ep_RA))**-1
example_delF = -np.log(ref_rep / example_reps)
F_ref = -np.log(pact) - np.log(ref_rep/Nns) + ep_RA;
deltaF = -np.log(ref_rep / rep_range) # in kT


# Load the actual data. 
fc_data = pd.read_csv('../../data/analyzed_foldchange.csv')
inferred_F = pd.read_csv('../../data/inferred_empirical_F.csv')

# Isolate the fc data to the relevant measurements 
fc_data = fc_data[(fc_data['temp']==37) &
                  (fc_data['strain']=='dilution') & 
                  (fc_data['repressors'] >= 0)].copy()
inferred_F = inferred_F[inferred_F['temp']==37].copy()

# Compute the summary statistics
rep_summary = fc_data.groupby(['date', 'run_number', 
                               'atc_ngml', 'carbon']).mean().reset_index()
summary = rep_summary.groupby(['atc_ngml', 'carbon']).agg(('mean', 'sem')).reset_index()

#%% Set up the figure axis
fig, ax = plt.subplots(2, 4, figsize=(7, 3.5), dpi=100)

# Define the axes and colors
ax_map = {'glucose':1, 'glycerol':2, 'acetate':3}
edgecolors = {'glucose':colors['dark_purple'], 'glycerol':colors['dark_green'],
              'acetate':colors['dark_brown']}
facecolors = {'glucose':colors['light_purple'], 'glycerol':colors['light_green'],
              'acetate':colors['light_brown']}

for g, d in summary.groupby(['carbon']):
    ax[0, ax_map[g]].errorbar(d['repressors']['mean'], d['fold_change']['mean'], 
                      xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'],
                      fmt='.', ms=8, color=edgecolors[g], 
                      markerfacecolor=facecolors[g], markeredgewidth=0.75)

for g, d in inferred_F.groupby(['carbon']):
        _reps = summary[(summary['carbon']==g)]
        F = d[d['parameter']=='empirical_F']['median']
        F_max = d[d['parameter']=='empirical_F']['hpd_max']
        F_min = d[d['parameter']=='empirical_F']['hpd_min']
        delF = F - F_ref   
        delF_max = F_max - F_ref   
        delF_min = F_min - F_ref   
        ax[1, ax_map[g]].errorbar(_reps['repressors']['mean'], -delF, xerr=_reps['repressors']['sem'], 
        fmt='.', color=edgecolors[g], markerfacecolor=facecolors[g], linestyle='none', linewidth=1,
        ms=8, markeredgewidth=0.75)
        ax[1, ax_map[g]].vlines(_reps['repressors']['mean'], -delF_max, -delF_min, lw=0.75,
                        color=edgecolors[g])

for i in range(4):
    ax[0, i].plot(rep_range, fc, 'k-')
    ax[0, i].set_xscale('log')
    ax[1, i].set_xscale('log')
    ax[1, i].plot(rep_range, deltaF, 'k-')
    ax[0, i].set_yscale('log')
    ax[0, i].set_ylabel('fold-change', style='italic', fontsize=8)
    ax[1, i].set_ylabel('∆F [k$_B$T]', style='italic', fontsize=8)

for a in ax.ravel():
    a.set_xlabel('repressors per cell', fontsize=8, style='italic')


# Plot the explanatory points. 
ax[0, 0].grid(False)
ax[1, 0].grid(False)
ax[0, 0].fill_between(rep_range, example_fc[2], 1.05, color=colors['light_red'],
            alpha=0.5)
ax[0, 0].fill_between(rep_range, example_fc[2], -.05, color=colors['light_green'],
            alpha=0.5)
ax[1, 0].fill_betweenx(np.linspace(-5, 5, 100), 100, 8000, 
                    color=colors['light_green'], alpha=0.5)
ax[1, 0].fill_betweenx(np.linspace(-5, 5, 100), 0, 99, 
                    color=colors['light_red'], alpha=0.5)
for i, r in enumerate(example_reps):
    ax[0, 0].plot(r, example_fc[i], 'o', markerfacecolor=example_facecolors[r],
            markeredgecolor=example_edgecolors[r], markeredgewidth=0.75)
    ax[1, 0].plot(r, example_delF[i], 'o', markerfacecolor=example_facecolors[r],
            markeredgecolor=example_edgecolors[r], markeredgewidth=0.75)

# ax[0, 0].grid(False)
plt.tight_layout()
plt.savefig('../../figs/Fig_delF_carbon_quality.pdf', bbox_inches='tight', 
            facecolor='white')
#%%
fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=100)
ax.set_xscale('log')
ax.set_yscale('log')

# Define the axes and colors
ax_map = {'glucose':1, 'glycerol':2, 'acetate':3}
edgecolors = {'glucose':colors['dark_purple'], 'glycerol':colors['dark_green'],
              'acetate':colors['dark_brown']}
facecolors = {'glucose':colors['light_purple'], 'glycerol':colors['light_green'],
              'acetate':colors['light_brown']}

ax.plot(rep_range, fc, 'k-')
for g, d in summary.groupby(['carbon']):
    ax.errorbar(d['repressors']['mean'], d['fold_change']['mean'], 
                      xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'],
                      fmt='.', ms=8, color=edgecolors[g], 
                      markerfacecolor=facecolors[g], markeredgewidth=0.25)
plt.savefig('../../figs/all_carbon_titration.pdf', bbox_inches='tight',
facecolor='white')
#%%
