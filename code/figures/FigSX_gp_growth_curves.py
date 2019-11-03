#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
colors, color_list = mwc.viz.personal_style()
bbox = dict(facecolor='none', edgecolor=colors['black'], lw=0.5)

# %%
plates = pd.read_csv('../../data/compiled_growth_plates.csv')
stats = pd.read_csv('../../data/compiled_growth_statistics.csv')
plates['temp_C'] = np.round(plates['temp_C'])
plates['time_min'] = np.round(plates['time_min'])
# plates = plates[plates['time_min'] <= 600]

# Load specific examples of GP processing. 
gp_output = pd.read_csv('../processing/growth_curves/20190206_r1_37C_mixedmedia_O2_growth/output/delta_glycerol/gp_output.csv')
plate = pd.read_csv('../processing/growth_curves/20190206_r1_37C_mixedmedia_O2_growth/output/growth_plate.csv')

# %%
# Instantiate figure
fig, ax = plt.subplots(1, 2, figsize=(5, 2.5))

# Add axis labels
ax[0].set_xlabel('time [min]', fontsize=8)
ax[0].set_ylabel('optical density (600 nm)', fontsize=8)
ax[1].set_xlabel('time [min]', fontsize=8)
ax[1].set_ylabel('growth rate [hr$^{-1}$]')

# Plot the growth curves and growth rates
ax[0].plot(gp_output['time'] - gp_output['time'].min(), gp_output['OD_raw_data'], '.', 
        markerfacecolor=colors['light_green'], 
        markeredgecolor=colors['green'], alpha=0.25, markeredgewidth=0.5, 
        label='raw measurements', ms=3)
ax[0].plot(gp_output['time'] - gp_output['time'].min(), 
        np.exp(gp_output['log(OD)_fit']), '-', color=colors['black'], lw=0.75,
        label='estimated value')
ax[1].fill_between(gp_output['time'] - gp_output['time'].min(),
            (gp_output['growth_rate'] + gp_output['growth_rate_std']) * 60, 
            (gp_output['growth_rate'] - gp_output['growth_rate_std']) * 60,
            color=colors['light_purple'], label='standard deviation', alpha=0.4)
ax[1].plot(gp_output['time'] - gp_output['time'].min(),
            gp_output['growth_rate'] * 60, '-', color=colors['dark_purple'], 
            label='mean')
# Add the legends
ax[1].legend()       
ax[0].legend()
mwc.viz.titlebox(ax[0], 'growth curve', size=8, color=colors['black'])
mwc.viz.titlebox(ax[1], 'growth rate', size=8, color=colors['black'])
plt.tight_layout()

# Add panel labels
fig.text(0.009, 0.9, '(A)', fontsize=10)
fig.text(0.5, 0.9, '(B)', fontsize=10)

# Add titles
plt.savefig('../../figs/FigSX_gp_growth_fig.pdf', facecolor='white', dpi=300)