import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pymc3 as pm
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


#%% Noise-free simulation.
# Set up parameters
num_div = 100
alpha_true = 150

# Set the protein number range.
prot_range = np.arange(10, 1000, 5)

# Set up a dataframe.
dfs = []
cols = ['summed', 'fluct', 'n1', 'n2', 'ntot']
for i, ntot in enumerate(prot_range):
    n1 = np.random.binomial(ntot, p=0.5, size=num_div)
    n2 = ntot - n1
    summed = alpha_true * ntot
    fluct = (alpha_true * n1 - alpha_true * n2)**2
    copy_dict = dict(summed=summed, fluct=fluct, n1=n1, n2=n2, ntot=ntot)
    df = pd.DataFrame(np.array([np.ones(num_div) * summed, fluct, n1, n2,
                                np.ones(num_div) * ntot]).T, columns=cols)
    dfs.append(df)
df = pd.concat(dfs, ignore_index=True)

#%% Plot the simulation data.
# Compute the averages.
grouped = df.groupby('ntot').mean()

# Compute the theory
theory = alpha_true * (alpha_true * prot_range)
fig, ax = plt.subplots(figsize=(5, 4))
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('$I_1 + I_2$')
ax.set_ylabel('$(I_1 + I_2)^2$')
_ = ax.plot(df['summed'], df['fluct'], '.', ms=0.5, label='simulated data')
_ = ax.plot(grouped['summed'], grouped['fluct'], '.', ms=4, label='averages')
_ = ax.plot(alpha_true * prot_range, theory, label='prediction')
_ = ax.legend()
sns.despine(offset=7)
plt.tight_layout()
plt.savefig('simulated_dilution_simple.pdf', bbox_inches='tight')

stats
# %%
with pm.Model() as model:
    like = mwc.bayes.DeterminsticCalibrationFactor(
        'alpha', alpha_true * df['n1'].values),
                                                   alpha_true * \
                                                       df['n2'].values,
                                                   testval = 1)
    trace=pm.sample(draws = 10000, tune = 10000)
    trace_df=mwc.stats.trace_to_dataframe(trace, model)
    stats=mwc.stats.compute_statistics(trace_df)

log_post=np.abs(log_post)
x=log_post / log_post[-1]

post=np.exp(-x) / np.sum(np.exp(-x))

# %%
fig, ax=plt.subplots(2, 1, figsize = (4.7, 3))
ax[0].set_ylabel(r'$\propto g(\alpha\, \vert \, [I_1, I_2])$')
ax[0].set_xlabel(r'$\alpha$ [a.u. / molecule]')
ax[1].set_ylabel(r'$\alpha$ [a.u. / molecule]')
ax[1].set_xlabel('step number')

_=ax[0].hist(trace_df['alpha'], bins = 100, alpha = 0.8, edgecolor = 'k',
               histtype='stepfilled', linewidth=0.8, normed=True,
               label='__nolegend__')
ypos = ax[0].get_ylim()[1] / 2
_ = ax[0].plot(stats['mode'], ypos, 'o', color='firebrick', ms=5, label='mode')
_ = ax[0].plot((stats['hpd_min'], stats['hpd_max']),(ypos, ypos),
               color='firebrick', ms=5)
_ = ax[1].plot(np.arange(0, len(trace_df), 1), trace_df['alpha'], alpha=0.6)
ax[0].text(-0.18, 1.2, '(A)', fontsize=8, transform=ax[0].transAxes)
ax[1].text(-0.18, 1.2, '(B)', fontsize=8, transform=ax[1].transAxes)
sns.despine(offset=7)
plt.tight_layout()
plt.savefig('alpha_simple_mcmc.pdf', bbox_inches='tight')

# %% Model 1
sigma = 0.001 * alpha_true

# Set the protein number range.
prot_range = np.arange(10, 1000, 5)

# Set up a dataframe.
dfs = []
cols = ['summed', 'fluct', 'I_1', 'I_2', 'n1', 'n2', 'ntot']
for i, ntot in enumerate(prot_range):
    n1 = np.random.binomial(ntot, p=0.5, size=num_div)
    n2 = ntot - n1
    I_1 = alpha_true * n1 + np.random.normal(0, sigma, size=num_div)
    I_2 = alpha_true * n2 + np.random.normal(0, sigma, size=num_div)
    summed = I_1 + I_2
    fluct = (I_1 - I_2)**2
    df = pd.DataFrame(np.array([summed, fluct, I_1, I_2, n1, n2,
                                np.ones(num_div) * ntot]).T, columns=cols)
    dfs.append(df)
model1_df = pd.concat(dfs, ignore_index=True)


# Perform the parameter estimation.
with pm.Model() as model1:
    alpha = pm.Uniform('alpha', lower=0, upper=1000)
    sigma = mwc.bayes.Jeffreys('sigma', lower=0.001 * alpha_true, upper=0.1 * alpha_true)
    n1 = pm.Uniform('n1', lower=1, upper=2000, shape=len(model1_df))
    ntot = pm.Uniform('ntot', lower=1, upper=2000, shape=len(model1_df))


    # Compute the means.
    mu1 = alpha * n1
    mu2 = alpha * (ntot - n1)

    # Define the likelihoods.
    like_1 = pm.Normal('like1', mu1, sd=sigma, observed=model1_df['I_1'].values)
    like_2 = pm.Normal('like2', mu2, sd=sigma, observed=model1_df['I_2'].values)

    # Burn with metropolis.
    step = pm.Metropolis()
    burn = pm.sample(tune=10000, draws=10000, step=step)

    # Perform the sampling.
    step = pm.NUTS()
    model1_trace = pm.sample(tune=1000, draws=1000, step=step, start=burn[-1])
    model1_trace_df = mwc.stats.trace_to_dataframe(model1_trace, model1)
    model1_stats = mwc.stats.compute_statistics(model1_trace_df)

# %%
