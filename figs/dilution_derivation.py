import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pymc3
import scipy.stats
import scipy.special
import sys
sys.path.insert(0, '../')
import mwc.bayes
import mwc.viz
import imp
imp.reload(mwc.viz)
mwc.viz.personal_style()
%matplotlib inline


# Set up parameters
num_div = 1000
alpha_seed = 150

# Distribute proteins via gamma.
n_tot = np.random.gamma(2, 40, size=num_div).astype(int)

# Binomially partition into daughters.
n_1 = np.random.binomial(n_tot, p=0.5)
n_2 = n_tot - n_1

# Convert to intensity using the alpha seed value.
I_1 = alpha_seed * n_1
I_2 = alpha_seed * n_2

# Assemble data in a dataframe for easy manipulation.
simple_df = pd.DataFrame(np.array((I_1, I_2, I_1 + I_2, (I_1 - I_2)**2)).T,
                         columns=['I_1', 'I_2', 'summed', 'fluct'])
simple_df.insert(3, 'alpha', alpha_seed)


# Bin the data by 50 events.
sorted_vals = simple_df.sort_values('summed')
bin_size = 50
bins = np.arange(0, len(sorted_vals) + bin_size, bin_size)

# Compute the means.
mean_summed, mean_fluct = [], []
for i in range(1, len(bins)):
    d = sorted_vals.iloc[bins[i - 1]:bins[i] + 1].mean()
    mean_summed.append(d['summed'])
    mean_fluct.append(d['fluct'])

# Compute the theoretical prediction.
I_tot_range = np.logspace(2, 5, 500)

theory = alpha_seed * I_tot_range

# Set up the figure canvas.
fig, ax = plt.subplots(1, 1, figsize=(4, 3))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('$I_1 + I_2$ [a.u.]')
ax.set_ylabel('$(I_1 - I_2)^2$ [a.u.]')

# Plot the data, bins, and theory.
_ = ax.plot(simple_df['summed'], simple_df['fluct'], '.', ms=1,
            label='simulated data')
_ = ax.plot(mean_summed, mean_fluct,  '.', ms=5, label='binned data')
_ = ax.plot(I_tot_range, theory, '-', label=r'$\alpha I_\mathrm{tot}$')

# Format
ax.set_xlim([500, 1E5])
_ = ax.legend()
# sns.despine(offset=5)
