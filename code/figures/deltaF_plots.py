#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import bokeh.io
import bokeh.plotting
import mwc.viz
import mwc.model
colors, color_list = mwc.viz.bokeh_theme()
mwc.viz.personal_style()


# Load the analyzed fold-change data
data = pd.read_csv('../../data/analyzed_foldchange.csv')

# Filter the data (drop )
data = data[(data['strain']=='dilution') & (data['repressors'] > 10) & (data['fold_change'] >= 0)]


#%%
# Define the reference state 
ref_rep = 100 # ARrbitrary choice
EP_RA = -13.9 # in kT
IPTG = 0 # in µM
EP_AI = 4.5 # in kT
TEMP_REF = 37 + 273 # in K
bohr_ref = mwc.model.SimpleRepression(R=ref_rep, ep_r=EP_RA, ep_ai=EP_AI,
                        effector_conc=IPTG, ka=139, ki=0.53).bohr_parameter()

# Compute what the free energy should be for each cell in the fold-change data
data['theo_deltaF'] = -np.log(data['repressors'] / ref_rep) 

# Compute the empirical bohr
data['empirical_bohr'] = np.log((1 / data['fold_change'].values) - 1)
data['empirical_deltaF'] = data['empirical_bohr'] - bohr_ref

#%%
# Start with only the carbon variants
carb = data[data['temp']==37]
temp = data[data['carbon']=='glucose']
# Compute the theory curve
rep_range = np.logspace(0, 6, 500)
theo_bohr = -np.log(ref_rep / rep_range)

# Set up a figure canvas
fig, ax = plt.subplots(2, 3, figsize=(6, 4), dpi=150) 
for i, a in enumerate(ax.ravel()):
    a.set_xscale('log')
    a.set_xlim([10, 2000])
    a.set_ylim([-6, 6])
    a.set_xlabel('repressors per cell')
    a.set_ylabel('$\Delta F$ [$k_BT$]')

    if i <=2:
        a.plot(rep_range, -np.log(ref_rep / rep_range), 'k-', lw=0.5, zorder=100)
    
# Set a dummy legend
ax[0,0].plot([], [], 'k-',  color='#C2C2C2', label='theory')
ax[0,0].plot([], [], '.', ms=2, color='#C2C2C2', label='single cell')
ax[0,0].errorbar([], [], [], fmt='.', ms=2, markeredgecolor='k', markeredgewidth=1,
                        markerfacecolor='#C2C2C2', label='binned by ATC',
                        linestyle='none', capsize=1, lw=1)
# ax[0, 0].legend(ncol=3, bbox_to_anchor=(0.5, 0.5))

# Assign the axes
carb_ax = {'glucose':ax[0,0], 'glycerol':ax[0, 1], 'acetate':ax[0, 2]}
temp_ax = {42: ax[1, 0], 37:ax[1, 1], 32:ax[1, 2]}
carb_colors = {'glucose':colors['purple'], 'glycerol':colors['green'],
                'acetate':colors['orange']}
temp_colors = {37: colors['purple'], 32: colors['blue'], 42: colors['red']}
for g, d in carb.groupby(['carbon']):

    d = d.copy()
    _ax = carb_ax[g]
    c = carb_colors[g]
    _ax.plot(d['repressors'], d['empirical_deltaF'], '.', ms=0.3, 
                rasterized=True, color='#C2C2C2', label='single cell data')

    # Round and group by repressor number
    grouped = d.groupby(['atc_ngml', 'date']).mean().reset_index()
    grouped = grouped.groupby(['atc_ngml']).agg(('mean', 'sem')).reset_index()
    _ax.errorbar(grouped['repressors']['mean'], grouped['empirical_deltaF']['mean'],
                xerr=grouped['repressors']['sem'], yerr=grouped['empirical_deltaF']['sem'], fmt='o', ms=3, 
                color=c, linestyle='none', lw=0.75, capsize=1, label='induction condition',
                markeredgecolor='k', markeredgewidth=0.15, zorder=99)
    _ax.set_title(f'{g}, 37°C', loc='left', style='italic')




for g, d in temp.groupby(['temp']):
    _ax = temp_ax[g]
    T_exp = g + 273
    delta_T = TEMP_REF / T_exp
    theo = -np.log(ref_rep/rep_range) + EP_RA * (1 - (TEMP_REF/T_exp)) - np.log((1 + np.exp(-EP_AI * delta_T)) / (1 + np.exp(-EP_AI)))
    _ax.plot(rep_range, theo, 'k-', lw=0.5, zorder=100)
    d = d.copy()
    c = temp_colors[g]
    _ax.plot(d['repressors'], d['empirical_deltaF'], '.', ms=0.3, 
                rasterized=True, color='#C2C2C2', label='single cell data')

    # Round and group by repressor number
    grouped = d.groupby(['atc_ngml', 'date', 'run_number']).mean().reset_index()
    grouped = grouped.groupby(['atc_ngml']).agg(('mean', 'sem')).reset_index()
    _ax.errorbar(grouped['repressors']['mean'], grouped['empirical_deltaF']['mean'],
                xerr=grouped['repressors']['sem'], yerr=grouped['empirical_deltaF']['sem'], fmt='o', ms=3, 
                color=c, linestyle='none', lw=0.75, capsize=1, label='induction condition',
                markeredgecolor='k', markeredgewidth=0.25, zorder=99)
    _ax.set_title(f'glucose, {g}°C', loc='left', style='italic')

plt.tight_layout()
plt.savefig('deltaF_variations.pdf', facecolor='w', bbox_inches='tight')
#%% Make the fold-change plots
rep_range = np.logspace(0, 4)
ref_fc = mwc.model.SimpleRepression(rep_range, ep_r=EP_RA, ep_ai=4.5, ka=139, ki=0.53,
                                    effector_conc=0).fold_change()

fig, ax = plt.subplots(2, 3, figsize=(6, 4), dpi=150)
for i, a in enumerate(ax.ravel()):
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlim([8, 2000])
    a.set_ylim([1E-4, 2])
    a.set_xlabel('repressors per cell')
    a.set_ylabel('fold-change')
    a.plot(rep_range, ref_fc, 'k-', lw=1, zorder=100)

carb_ax = {'glucose':ax[0,0], 'glycerol':ax[0, 1], 'acetate':ax[0, 2]}
temp_ax = {42: ax[1, 0], 37:ax[1, 1], 32:ax[1, 2]}
carb_colors = {'glucose':colors['purple'], 'glycerol':colors['green'],
                'acetate':colors['orange']}
temp_colors = {37: colors['purple'], 32: colors['blue'], 42: colors['red']}

for g, d in carb.groupby(['carbon']):

    d = d.copy()
    _ax = carb_ax[g]
    c = carb_colors[g]
    _ax.plot(d['repressors'], d['fold_change'], '.', ms=0.3, 
                rasterized=True, color='#C2C2C2', label='single cell data')

    # Round and group by repressor number
    grouped = d.groupby(['atc_ngml', 'date']).mean().reset_index()
    grouped = grouped.groupby(['atc_ngml']).agg(('mean', 'sem')).reset_index()
    _ax.errorbar(grouped['repressors']['mean'], grouped['fold_change']['mean'],
                xerr=grouped['repressors']['sem'], yerr=grouped['fold_change']['sem'], fmt='o', ms=3, 
                color=c, linestyle='none', lw=0.75, capsize=1, label='induction condition',
                markeredgecolor='k', markeredgewidth=0.15, zorder=99)
    _ax.set_title(f'{g}, 37°C', loc='left', style='italic')



for g, d in temp.groupby(['temp']):
    _ax = temp_ax[g]
    T_exp = g + 273
    delta_T = TEMP_REF / T_exp
    theo = mwc.model.SimpleRepression(rep_range, ep_r=EP_RA * delta_T, 
                        ka=139, ki=0.53, effector_conc=0, ep_ai=EP_AI * delta_T).fold_change()
    _ax.plot(rep_range, theo, temp_colors[g], lw=0.5, zorder=100)
    d = d.copy()
    c = temp_colors[g]
    _ax.plot(d['repressors'], d['fold_change'], '.', ms=0.3, 
                rasterized=True, color='#C2C2C2', label='single cell data')

    # Round and group by repressor number
    grouped = d.groupby(['atc_ngml', 'date', 'run_number']).mean().reset_index()
    grouped = grouped.groupby(['atc_ngml']).agg(('mean', 'sem')).reset_index()
    _ax.errorbar(grouped['repressors']['mean'], grouped['fold_change']['mean'],
                xerr=grouped['repressors']['sem'], yerr=grouped['fold_change']['sem'], fmt='o', ms=3, 
                color=c, linestyle='none', lw=0.75, capsize=1, label='induction condition',
                markeredgecolor='k', markeredgewidth=0.25, zorder=99)
    _ax.set_title(f'glucose, {g}°C', loc='left', style='italic')


plt.tight_layout() 
plt.savefig('foldchange_variations.pdf', bbox_inches='tight', facecolor='w')
#%%
