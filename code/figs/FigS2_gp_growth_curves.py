"""
Author:
        Griffin Chure
License: 
        MIT
Description:
        Generates a figure with three subplots showing the raw data from a
        growth curve, a fit from Gaussian process modeling, and the plot of
        maximum growth rate vs time.
Required Data Sets:
        glycerol_example_gp_output.csv
        example_growth_plate.csv
"""
#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
colors, color_list = mwc.viz.personal_style()

# Load specific examples of GP processing. 
gp_output = pd.read_csv('../data/glycerol_example_gp_output.csv', comment='#')
plate = pd.read_csv('../data/example_growth_plate.csv', comment='#')
glyc = plate[(plate['carbon']=='glycerol') & (plate['time_min']<1200)]

# %%
# Instantiate figure
fig, ax = plt.subplots(1, 3, figsize=(7, 2.5))

# Add axis labels
for i in range(2):
        ax[i].set_xlabel('time [min]', fontsize=8)
        ax[i].set_ylabel('optical density (600 nm)', fontsize=8)
ax[0].set_ylim([0, 0.45])
ax[2].set_xlabel('time [min]', fontsize=8)
ax[2].set_ylabel('growth rate [hr$^{-1}$]')

# Plot the growth curves and growth rates
ax[0].plot(glyc['time_min'], glyc['od_sub'], '.', color=colors['green'],
           markeredgewidth=0.5, ms=1, alpha=0.15, label='__nolegend__')
ax[0].fill_betweenx(np.linspace(0, 0.6, 100), gp_output['time'].min(), gp_output['time'].max(), 
                    color=colors['light_grey'], alpha=0.25, label='exponential region')
ax[1].plot(gp_output['time'], gp_output['OD_raw_data'], '.', 
        markerfacecolor=colors['light_green'], 
        markeredgecolor=colors['green'], alpha=0.5, markeredgewidth=0.1, 
        label='__nolegend__', ms=3)
ax[1].plot(gp_output['time'], 
        np.exp(gp_output['log(OD)_fit']), '-', color=colors['purple'], lw=0.75,
        label='estimated value')
ax[2].fill_between(gp_output['time'],
            (gp_output['growth_rate']) * 60, 
            (gp_output['growth_rate']) * 60,
            color=colors['light_purple'], label='standard deviation', alpha=0.4)
ax[2].plot(gp_output['time'],
            gp_output['growth_rate'] * 60, '-', color=colors['dark_purple'], 
            label='mean')
# Add the legends
ax[0].legend()
ax[1].legend()       
ax[2].legend()

# Add the titleboxes
mwc.viz.titlebox(ax[0], 'GROWTH CURVE', size=6, boxsize="10%", color=colors['black'])
mwc.viz.titlebox(ax[1], 'EXPONENTIAL REGION', size=6, boxsize="10%", color=colors['black'])
mwc.viz.titlebox(ax[2], 'GROWTH RATE', size=6, boxsize="10%", color=colors['black'])
plt.tight_layout()

# Add panel labels
fig.text(0.009, 0.9, '(A)', fontsize=9, fontweight='bold')
fig.text(0.35, 0.9, '(B)', fontsize=9, fontweight='bold')
fig.text(0.67, 0.9, '(C)', fontsize=9, fontweight='bold')

# Add titles
plt.savefig('../../figs/FigS2_gp_growth_fig.pdf', facecolor='white', dpi=300)

# %%
