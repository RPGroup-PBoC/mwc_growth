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
import tqdm
import imp
imp.reload(mwc.viz)
colors = mwc.viz.personal_style()
%matplotlib inline
np.random.seed(666)

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
# Estimate parameter through minimization, compute posterior and approximation.
est, err = mwc.bayes.estimate_calibration_factor(
    df['n1'] * alpha_true, df['n2'] * alpha_true)
alpha_range = np.linspace(est - 20, est + 20, 500)
log_post = np.zeros_like(alpha_range)
for i, a in enumerate(alpha_range):
    log_post[i] = mwc.bayes.deterministic_log_posterior(a, df['n1'] * alpha_true,
                                                        df['n2'] * alpha_true)

# Normalize the log posterior and compute.
posterior = np.exp(log_post - scipy.special.logsumexp(log_post))
approx = scipy.stats.norm.pdf(alpha_range, loc=est, scale=err)
approx = approx / np.sum(approx)
# %% Plot the posterior and approximation.
fig, ax = plt.subplots(1, 1, figsize=(5, 4))
ax.set_xlabel('calibration factor [a.u. / molecule]')
ax.set_ylabel('probability')

_ = ax.plot(alpha_range, posterior, label='posterior')
_ = ax.fill_between(alpha_range, posterior, alpha=0.5, label='__nolegend__')
_ = ax.plot(alpha_range, approx, ':', label='Gaussian approximation', lw=1.5)
_ = ax.vlines(alpha_true, np.min(approx),
              np.max(approx), label=r'true $\alpha$', color=colors[4],
              zorder=1000)
_ = ax.legend()
mwc.viz.format_axes()
plt.tight_layout()
plt.savefig('alpha_simple_minimization.pdf', bbox_inches='tight')

est, err
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

# %% Measured alpha as a function of noise.
%matplotlib inline
noise_range = np.logspace(-3, 2, 500) * alpha_true
alpha_bar = np.empty_like(noise_range)
for j, sigma in enumerate(tqdm.tqdm(noise_range)):
    _df = []
    for n in prot_range:
        n1 = np.random.binomial(n, 0.5, size=num_div)
        n2 = n - n1

        I_1 = alpha_true * n1 + \
            np.random.normal(0, sigma, size=num_div)
        I_2 = alpha_true * n2 + \
            np.random.normal(0, sigma, size=num_div)
        __df = pd.DataFrame(np.array([I_1, I_2]).T, columns=['I_1', 'I_2'])
        _df.append(__df)
    _df = pd.concat(_df, ignore_index=True)
    _df = _df[_df >= 0]

    # Estimate alpha bar via minimization.
    x, s = mwc.bayes.estimate_calibration_factor(_df['I_1'], _df['I_2'])
    alpha_bar[j] = x
# %%
fig, ax = plt.subplots(1, 1, figsize=(4.5, 3))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\sigma_1 / \alpha_{true}$')
ax.set_ylabel(r'$\alpha_{est} / \alpha_{true}$')
ax.set_title('random measurement error', fontsize=10)
_ = ax.plot(noise_range / alpha_true, alpha_bar / alpha_true,
            '.', ms=1, label='estimated\nvalue')

# # Plot the approximate errors for illumination.
max_height = np.max(alpha_bar / alpha_true)
_ = ax.plot([], [], lw=2, alpha=0.5,
            label='Hg lamp', color='dodgerblue')
_ = ax.plot([], [], lw=2, alpha=0.5, label='LED', color=colors[1])
_ = ax.plot([], [], lw=2, alpha=0.5,
            label='laser', color=colors[3])
_ = ax.legend(loc='upper right')

# Make the inset.
mwc.viz.format_axes()
ax2 = plt.axes([0.2, 0.5, 0.3, 0.3], facecolor='w')
ax2.plot(noise_range / alpha_true, alpha_bar / alpha_true, '.',
         ms=1)
ax2.set_xlim([0.01, 0.025])
ax2.set_ylim([0.95, 1.05])

max_height = np.max(alpha_bar / alpha_true)
_ = ax2.vlines(0.014, 0.95, 1.05, lw=3, alpha=0.5,
               label='Hg lamp', color='dodgerblue')
_ = ax2.vlines(0.014, 0.95, 1.05, lw=2, alpha=0.5,
               label='LED', color=colors[1])
_ = ax2.vlines(0.02, 0.95, 1.05, lw=3, alpha=0.5,
               label='laser', color=colors[3])
ax2.grid(False)
_ = ax.legend(loc='upper right')
plt.savefig('error_est_model1.pdf', bbox_inches='tight')
# plt.tight_layout()
# plt.savefig('error_est_model1.pdf', bbox_inches='tight')

# %% Model II.
noise_range = np.logspace(-4, 0, 500)
alpha_bar = np.empty_like(noise_range)
alpha_err = np.empty_like(noise_range)
for j, sigma in enumerate(noise_range):
    dfs = []
    for i, ntot in enumerate(prot_range):
        n1 = np.random.binomial(ntot, p=0.5, size=num_div)
        n2 = ntot - n1
        phi = np.random.normal(1, sigma, size=num_div)
        I_1 = n1 * alpha_true * phi
        I_2 = n2 * alpha_true * phi
        _df = pd.DataFrame(np.array([I_1, I_2]).T, columns=['I_1', 'I_2'])
        dfs.append(_df)
    _df = pd.concat(dfs, ignore_index=True)

    # Drop negative values for proper estimation.
    _df = _df[_df > 0]

    # Estimate alpha bar via minimization.
    x, s = mwc.bayes.estimate_calibration_factor(_df['I_1'], _df['I_2'])
    alpha_bar[j] = x

# %%
fig, ax = plt.subplots(1, 1, figsize=(4.5, 3))
ax.set_xscale('log')
ax.set_xlabel('$\sigma_2$')
ax.set_ylabel(r'$\alpha_{est} / \alpha_{true}$')
_ = ax.plot(noise_range, alpha_bar / alpha_true,
            '.', ms=2, label='estimated\nvalue')
mwc.viz.format_axes()
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
