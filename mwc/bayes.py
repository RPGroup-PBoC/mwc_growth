# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pymc3 as pm
import scipy.special
import scipy.optimize
import statsmodels.tools.numdiff as smnd
import theano.tensor as tt


def log_prior():
    """
    Computes the log prior for alpha and I_2 as needed in the model. These are
    taken to be uniform on the range of the camera bitdepth, 0 - 65535.
    """
    return -32 * np.log(2)


def deterministic_log_posterior(alpha, I_1, I_2, p=0.5, neg=False):
    """
    Computes the log posterior of the deterministic model for the calibration
    factor.

    Parameters
    ----------
    alpha : float
        The calibration factor in units of a.u. per molecule. This must be
        positive
    I_1, I_2 : 1d-arrays or Pandas Series.
        The intensity of the two sister cells in units of a.u. per cell.
        Negative values will raise a ValueError.
    p: float between 0 and 1
        The partitioning probability into one cell or the other. Default value
        is fair partitioning, 0.5.
    neg : bool
        If True, the negative log posterior is returned. Default is False

    Returns
    -------
    logp : float
        Value of the log posterior with the provided parameter values.
    """
    # Determine the prefactor.
    if neg is True:
        prefactor = -1
    else:
        prefactor = 1

    # Ensure alpha is positive. If not, return
    if alpha < 0:
        return prefactor * -np.inf

    # Ensure that the two intensities are positive.
    if (I_1 < 0).any() or (I_2 < 0).any():
        raise ValueError('I_1 or I_2 contains negative values. Fix that plz.')

    # Make sure value for p is sensical.
    if (p < 0) | (p > 1):
        raise ValueError('p must be on the domain [0, 1]')

    # Convert the intensities to protein number.
    n_1 = I_1 / alpha
    n_2 = I_2 / alpha
    n_tot = n_1 + n_2
    k = len(I_1)

    # Compute the various parts of the posterior.
    binom = scipy.special.gammaln(n_tot + 1).sum() - scipy.special.gammaln(
        n_1 + 1).sum() - scipy.special.gammaln(n_2 + 1).sum()
    prob = n_1.sum() * np.log(p) + n_2.sum() * np.log(1 - p)
    change_of_var = -k * np.log(alpha)

    # Compute the log prior.
    logprior = log_prior()

    # Assemble the log posterior.
    logpost = change_of_var + binom + prob
    return prefactor * (logpost + logprior)


def estimate_calibration_factor(I_1, I_2, p=0.5, return_eval=False):
    """
    Estimates the fluorescence calibration factor through optimization.

    Parameters
    ----------
    I_1, I_2 : 1d-arrays or Pandas Series
        The intensities of two sister cells.
    p : float
        The probability of partitioning into one cell or another. Default is
        fair partitioning (0.5).
    return_eval : Bool
        If True, the evaluation statistics from the optimization will be returned.

    Returns
    -------
    alpha_opt, alpha_std : float
        The best-fit value for alpha and standard devation.
    """

    # Perform data validation checks.
    if (I_1 < 0).any() | (I_2 < 0).any():
        raise ValueError(
            'I_1 and I_2 may not contain negative values. Fix that plz.')
    if (p < 0) | (p > 1):
        raise ValueError('p must be between 0 and 1.')

    # Perform the optimization
    popt = scipy.optimize.minimize_scalar(
        deterministic_log_posterior, args=(I_1, I_2, p, True))
    alpha_opt = popt.x

    # Compute the hessian.
    hess = smnd.approx_hess([alpha_opt], deterministic_log_posterior,
                            args=(I_1, I_2, p, False))
    cov = -np.linalg.inv(hess)
    alpha_std = np.sqrt(cov[0])[0]
    if return_eval is True:
        return [alpha_opt, alpha_std, popt]

    else:
        return [alpha_opt, alpha_std]

def chains_to_dataframe(fit, varnames=None):
    """
    Converts the generated traces from MCMC sampling to a tidy
    pandas DataFrame.

    Parameters
    ----------
    fit : pystan sampling output
        The raw MCMC output.
    var_names : list of str
        Names of desired parameters. If `None`, all parameters will be
        returned.

    Returns
    -------
    df : pandas DataFrame
        Pandas DataFrame containing all samples from the MCMC.
    """

    data = fit.extract()
    keys = list(data.keys())
    if varnames == None:
        varnames = [k for k in keys if 'lp__' not in k]
    else:
        varnames = fit.unconstrained_param_names()

    samples = {}
    for i, key in enumerate(varnames):
        # Get the shape.
        dim = np.shape(data[key])
        if len(dim) == 2:
            for j in range(dim[-1]):
                samples['{}.{}'.format(key, j+1)] = data[key][:, j]
    
        else:
            samples[key] = data[key]
    samples['logp'] = data['lp__']
    return pd.DataFrame(samples)
    