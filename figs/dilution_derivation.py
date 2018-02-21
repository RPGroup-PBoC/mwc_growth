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
ax.set_ylabel('$(I_1 - I_2)^2$')
_ = ax.plot(df['summed'], df['fluct'], '.', ms=0.5, label='simulated data')
_ = ax.plot(grouped['summed'], grouped['fluct'], '.', ms=4, label='averages')
_ = ax.plot(alpha_true * prot_range, theory, label='prediction')
_ = ax.legend()
sns.despine(offset=7)
plt.tight_layout()
plt.savefig('simulated_dilution_simple.pdf', bbox_inches='tight')


# %%
with pm.Model() as model:
    like = mwc.bayes.DeterminsticCalibrationFactor(
        'alpha', alpha_true * df['n1'].values, alpha_true * df['n2'].values, testval=1)
    trace = pm.sample(draws=10000, tune=10000)
    trace_df = mwc.stats.trace_to_dataframe(trace, model)
    stats = mwc.stats.compute_statistics(trace_df)

# %%
fig, ax = plt.subplots(2, 1, figsize=(4.7, 3))
ax[0].set_ylabel(r'$\propto g(\alpha\, \vert \, [I_1, I_2])$')
ax[0].set_xlabel(r'$\alpha$ [a.u. / molecule]')
ax[1].set_ylabel(r'$\alpha$ [a.u. / molecule]')
ax[1].set_xlabel('step number')

_ = ax[0].hist(trace_df['alpha'], bins=100, alpha=0.8, edgecolor='k',
               histtype='stepfilled', linewidth=0.8, normed=True,
               label='__nolegend__')
ypos = ax[0].get_ylim()[1] / 2
_ = ax[0].plot(stats['mode'], ypos, 'o',
               color='firebrick', ms=5, label='mode')
_ = ax[0].plot((stats['hpd_min'], stats['hpd_max']), (ypos, ypos),
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

# Compute the means.
mean_vals = model1_df.groupby('ntot').mean()


# %% Measured alpha as a function of noise.
%matplotlib inline
noise_range = np.logspace(-5, -1, 500)
alpha_bar = np.empty_like(noise_range)
alpha_err = np.empty_like(noise_range)
for j, sigma in enumerate(noise_range):
    dfs = []
    for i, ntot in enumerate(prot_range):
        n1 = np.random.binomial(ntot, p=0.5, size=num_div)
        n2 = ntot - n1
        I_1 = alpha_true * n1 + \
            np.random.normal(0, sigma * alpha_true * n1, size=num_div)
        I_2 = alpha_true * n2 + \
            np.random.normal(0, sigma * alpha_true * n2, size=num_div)
        summed = I_1 + I_2
        fluct = (I_1 - I_2)**2
        df = pd.DataFrame(np.array([I_1, I_2]).T, columns=['I_1', 'I_2'])
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)

    # Estimate alpha bar via minimization.
    x, s = mwc.bayes.estimate_calibration_factor(df['I_1'], df['I_2'])
    alpha_bar[j] = x
    alpha_err[j] = s


# %%
fig, ax = plt.subplots(2, 1, figsize=(5, 4.5))
ax[0].set_xscale('log')
ax[0].set_xlabel('$\sigma$ (% of $I_{tot}$)')
ax[0].set_ylabel(r'$\alpha$ [a.u. / mol.]')
ax[1].set_xlabel('$\sigma$ (% of  $I_{tot}$)')
ax[1].set_ylabel(r'$\alpha_{est.} / \alpha_{true}$')
ax[1].set_xscale('log')
ax[0].set_title('model I - measurement noise', fontsize=10)
ax[0].set_ylim([100, 600])
ax[0].text(-0.15, 1.1, '(A)', fontsize=8, transform=ax[0].transAxes)
ax[1].text(-0.15, 1.1, '(B)', fontsize=8, transform=ax[1].transAxes)
_ = ax[0].plot(noise_range * 100, alpha_bar, '.',
               ms=1, label='estimated value')
_ = ax[0].hlines(alpha_true, noise_range[0] * 100,
                 noise_range[-1] * 100, label='true value', color='dodgerblue')
_ = ax[0].legend()
_ = ax[1].plot(noise_range * 100, alpha_bar / alpha_true,
               '.', ms=1, label='estimated value')


# Plot the approximate errors for illumination.
_ = ax[1].vlines(0.1, 1, 3.6, lw=2, alpha=0.5,
                 label='Hg lamp', color='firebrick')
_ = ax[1].vlines(0.2, 1, 3.6, lw=2, alpha=0.5, label='LED', color='slategray')
_ = ax[1].vlines(0.4, 1, 3.6, lw=2, alpha=0.5,
                 label='laser', color='dodgerblue')
_ = ax[1].legend()
sns.despine(offset=7)
plt.tight_layout()
plt.savefig('error_est_model1.pdf', bbox_inches='tight')

# %% Model II.
noise_range = np.logspace(-5, 0, 500)
alpha_bar = np.empty_like(noise_range)
alpha_err = np.empty_like(noise_range)
for j, sigma in enumerate(noise_range):
    dfs = []
    for i, ntot in enumerate(prot_range):
        n1 = np.random.binomial(ntot, p=0.5, size=num_div)
        n2 = ntot - n1
        I_1 = np.array(
            [np.sum(np.random.normal(alpha_true, sigma * alpha_true, size=n)) for n in n1])
        I_2 = np.array(
            [np.sum(np.random.normal(alpha_true, sigma * alpha_true, size=n)) for n in n2])
        summed = I_1 + I_2
        fluct = (I_1 - I_2)**2
        df = pd.DataFrame(np.array([I_1, I_2]).T, columns=['I_1', 'I_2'])
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)

    # Drop negative values for proper estimation.
    df = df[df > 0]

    # Estimate alpha bar via minimization.
    x, s = mwc.bayes.estimate_calibration_factor(df['I_1'], df['I_2'])
    alpha_bar[j] = x
    alpha_err[j] = s

# %%
fig, ax = plt.subplots(2, 1, figsize=(5, 4.5))
ax[0].set_xlabel('$\sigma$ (% of $I_{tot}$)')
ax[0].set_xscale('log')
ax[0].set_ylabel(r'$\alpha$ [a.u. / mol.]')
ax[1].set_xlabel('$\sigma$ (% of  $I_{tot}$)')
ax[1].set_ylabel(r'$\alpha_{est.} / \alpha_{true}$')
ax[0].set_title('model II - spatial variation', fontsize=10)
ax[1].set_xscale('log')
ax[0].text(-0.15, 1.1, '(A)', fontsize=8, transform=ax[0].transAxes)
ax[1].text(-0.15, 1.1, '(B)', fontsize=8, transform=ax[1].transAxes)
_ = ax[0].plot(noise_range * 100, alpha_bar, '.',
               ms=1, label='estimated value')
_ = ax[0].hlines(alpha_true, noise_range[0] * 100,
                 noise_range[-1] * 100, label='true value', color='dodgerblue')
_ = ax[0].legend()
_ = ax[1].plot(noise_range * 100, alpha_bar / alpha_true,
               '.', ms=1, label='estimated value')

# Plot the approximate errors for illumination.
_ = ax[1].legend()
sns.despine(offset=7)
plt.tight_layout()
plt.savefig('err_est_model2.pdf', bbox_inches='tight')

# %% Model III
noise_range = np.logspace(-5, 0, 500)
alpha_bar = np.empty_like(noise_range)
alpha_err = np.empty_like(noise_range)
for j, sigma in enumerate(noise_range):
    dfs = []
    for i, ntot in enumerate(prot_range):
        n1 = np.random.binomial(ntot, p=0.5, size=num_div)
        n2 = ntot - n1
        I_1 = np.random.normal(1, sigma, size=num_div) * alpha_true * n1
        I_2 = np.random.normal(1, sigma, size=num_div) * alpha_true * n2
        summed = I_1 + I_2
        fluct = (I_1 - I_2)**2
        df = pd.DataFrame(np.array([I_1, I_2]).T, columns=['I_1', 'I_2'])
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)

    # Drop negative values for proper estimation.
    df = df[df > 0]

    # Estimate alpha bar via minimization.
    x, s = mwc.bayes.estimate_calibration_factor(df['I_1'], df['I_2'])
    alpha_bar[j] = x
    alpha_err[j] = s

# %% Plotting
fig, ax = plt.subplots(2, 2, figsize=(5, 4.5))
ax[0, 0].set_xlabel('$\sigma$ (% of $I_{tot}$)')
ax[0, 1].set_xlabel('$\sigma$ (% of $I_{tot}$)')
ax[0, 0].set_xscale('log')
ax[0, 1].set_xscale('log')
ax[1, 1].set_xscale('log')
ax[1, 0].set_xscale('log')
ax[1, 0].set_yscale('log')
ax[0, 0].set_ylabel(r'$\alpha$ [a.u. / mol.]')
ax[0, 1].set_ylabel(r'$\alpha$ [a.u. / mol.]')
ax[0, 1].set_xlabel('$\sigma$ (% variation)')
ax[1, 1].set_xlabel('$\sigma$ (% variation)')
ax[0, 0].set_xlabel('$\sigma$ (% variation)')
ax[1, 0].set_xlabel('$\sigma$ (% variation)')
ax[1, 0].set_ylabel(r'$\alpha_{est.} / \alpha_{true}$')
ax[1, 1].set_ylabel(r'$\alpha_{est.} / \alpha_{true}$')
fig.text(0.4, 1., 'model III - temporal variation', fontsize=10)
# ax[0].text(-0.15, 1.1, '(A)', fontsize=8, transform=ax[0].transAxes)
# ax[1].text(-0.15, 1.1, '(B)', fontsize=8, transform=ax[1].transAxes)
_ = ax[0, 0].plot(noise_range * 100, alpha_bar, '.',
                  ms=1, label='estimated value')
_ = ax[0, 1].plot(noise_range * 100, alpha_bar, '.',
                  ms=1, label='estimated value')
_ = ax[0, 0].hlines(alpha_true, noise_range[0] * 100,
                    noise_range[-1] * 100, label='true value', color='dodgerblue')
_ = ax[0, 1].hlines(alpha_true, noise_range[0] * 100,
                    noise_range[-1] * 100, label='true value', color='dodgerblue')
_ = ax[0, 0].legend()
_ = ax[1, 0].plot(noise_range * 100, alpha_bar / alpha_true,
                  '.', ms=1, label='estimated value')

_ = ax[1, 1].plot(noise_range * 100, alpha_bar / alpha_true,
                  '.', ms=1, label='estimated value')

ax[0, 1].set_ylim([100, 750])
ax[1, 1].set_ylim([0.8, 3])
ax[0, 0].text(-0.5, 1.05, '(A)', fontsize=8, transform=ax[0, 0].transAxes)
ax[0, 1].text(-0.5, 1.05, '(B)', fontsize=8, transform=ax[1, 0].transAxes)
ax[1, 0].text(-0.5, 1.05, '(C)', fontsize=8, transform=ax[1, 0].transAxes)
ax[1, 1].text(-0.5, 1.05, '(D)', fontsize=8, transform=ax[1, 1].transAxes)
# Plot the approximate errors for illumination.
# Plot the approximate errors for illumination.
# _ = ax[1].legend()
sns.despine(offset=7)
plt.tight_layout()
plt.savefig('err_est_model3.pdf', bbox_inches='tight')


# %% Perform the parameter estimation.
with pm.Model() as model1:
    alpha = pm.Uniform('alpha', lower=0, upper=1000)
    sigma = mwc.bayes.Jeffreys('sigma', lower=0.001 *
                               alpha_true, upper=0.1 * alpha_true)
    n1 = pm.Uniform('n1', lower=1, upper=2000, shape=len(model1_df))
    ntot = pm.Uniform('ntot', lower=1, upper=2000, shape=len(model1_df))

    # Compute the means.
    mu1 = alpha * n1
    mu2 = alpha * (ntot - n1)

    # Define the likelihoods.
    like_1 = pm.Normal('like1', mu1, sd=sigma,
                       observed=model1_df['I_1'].values)
    like_2 = pm.Normal('like2', mu2, sd=sigma,
                       observed=model1_df['I_2'].values)

    # Burn with metropolis.
    step = pm.Metropolis()
    burn = pm.sample(tune=10000, draws=10000, step=step)

    # Perform the sampling.
    step = pm.NUTS()
    model1_trace = pm.sample(tune=1000, draws=1000, step=step, start=burn[-1])
    model1_trace_df = mwc.stats.trace_to_dataframe(model1_trace, model1)
    model1_stats = mwc.stats.compute_statistics(model1_trace_df)

# %%
