#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import mwc.bayes
import scipy.stats
import scipy.special
colors, _ = mwc.viz.personal_style()

# %%
# Load data, restrict, and define new constants
fluct_data = pd.read_csv('../../data/analyzed_fluctuations.csv')
fluct_data = fluct_data[fluct_data['date']==20181002]
alpha_range = np.linspace(550, 950, 500)

#%%
model = mwc.bayes.StanModel('../stan/calibration_factor.stan')
fit, samples = model.sample(dict(N=len(fluct_data), I1=fluct_data['I_1'], 
                                I2=fluct_data['I_2']))
#%%
gauss_approx = scipy.stats.norm(np.mean(samples['alpha']), 
                               np.std(samples['alpha'])).pdf(alpha_range)
#%% 
# Instantiate the figure
fig, ax = plt.subplots(2, 2, figsize=(6, 5))

# Format the axes
# ax[0].set_xscale('log')
# ax[0].set_yscale('log')
ax[0, 0].set_xlabel('$\log_{10}$ cell intensity')
ax[0, 0].set_ylabel('number of observations')
ax[0, 1].set_xlabel('$I_1 / (I_1 + I_2)$')
ax[0, 1].set_ylabel('number of observations')
# ax[0, 1].set_ylim([])
ax[1, 0].set_xlabel('cell volume [fL]')
ax[1, 0].set_ylabel('fractional intensity')
ax[1, 1].set_xlabel('calibration factor [a.u. / LacI]')
ax[1, 1].set_ylabel('probability')
ax[1, 1].set_xlim([550, 950])
ax[1, 0].set_xlim(0, 3.5)
mwc.viz.titlebox(ax[0, 0], 'partitioned intensity', color=colors['black'])
mwc.viz.titlebox(ax[0, 1], 'fractional intensity', color=colors['black'])
mwc.viz.titlebox(ax[1, 0], 'partitioning probability', color=colors['black'])
mwc.viz.titlebox(ax[1, 1], r'$\alpha$ posterior distribution', color=colors['black'])
avg_partitioning  = np.mean(fluct_data['I_1'] / fluct_data['summed'])

# Plot the fluctuations
# ax[0].plot(fluct_data['summed'], fluct_data['fluct'], 'k.', ms=1, alpha=0.75)
ax[0,0].hist(np.log10(fluct_data['I_1']), bins=50, color=colors['light_purple'],
            edgecolor=colors['dark_purple'], label='$I_1$')
ax[0,0].hist(np.log10(fluct_data['I_2']), bins=50, color=colors['orange'],
            edgecolor=colors['dark_orange'], alpha=0.5, label='$I_2$')

ax[0, 1].hist(fluct_data['I_1'] / fluct_data['summed'], bins=50, color='grey', 
             edgecolor=colors['black'], alpha=0.5, lw=0.25)
ax[0, 1].vlines(0.5, 0, ax[0, 1].get_ylim()[1], color=colors['orange'], label='even partitioning\n(0.50)')
ax[0, 1].vlines(avg_partitioning, 0, ax[0, 1].get_ylim()[1], 
        color=colors['dark_purple'], label=f'average value\n({np.round(avg_partitioning, decimals=3)})')
        
ax[1, 0].plot(fluct_data['volume_1'], fluct_data['I_1'] / fluct_data['summed'],
              'k.', ms=0.75, alpha=0.75)
ax[1, 0].plot(fluct_data['volume_2'], fluct_data['I_2'] / fluct_data['summed'],
              'k.', ms=0.75,  alpha=0.75)

ax[1, 1].hist(samples['alpha'], bins=25, color=colors['purple'], edgecolor=colors['dark_purple'],
            lw=0.1, density=True, alpha=0.75, label='sampled\ndistribution')
ax[1, 1].plot(alpha_range, gauss_approx, color=colors['orange'], 
            label='gaussian\napproximation', lw=1.5)
ax[0, 0].legend()
ax[0, 1].legend(handlelength=1, fontsize=5.5)
ax[1, 1].legend(handlelength=1, fontsize=5.5)
plt.tight_layout()
fig.text(0, 0.95, '(A)', fontsize=10)
fig.text(0.5, 0.95, '(B)', fontsize=10)
fig.text(0, 0.46, '(C)', fontsize=10)
fig.text(0.5, 0.46, '(D)', fontsize=10)
plt.savefig('../../figs/FigSX_partitioning_statistics.pdf', bbox_inches='tight',
            facecolor='white')
# %%


# %%
