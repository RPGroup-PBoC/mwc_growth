# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pymc3 as pm
import scipy.special
import scipy.optimize
import statsmodels.tools.numdiff as smnd
import theano.tensor as tt


class MarginalizedHomoscedasticNormal(pm.Continuous):
    """
    A bivariate Normal distribution after marginalization of the
    variance, sigma.

    .. math::

        g(\mu, k \vert y) = \left((y - \mu)^2\right)^{-k/2}


    Parameters
    ----------
    mu : PyMC3 RV object
        Mean components of the distribution.
    """

    def __init__(self, mu=None, *args, **kwargs):
        super(MarginalizedHomoscedasticNormal, self).__init__(*args, **kwargs)
        self.mu = mu = pm.theanof.floatX(tt.as_tensor_variable(mu))
        self.median = mu
        self.mode = mu
        self.mean = mu

    def logp(self, values):
        k = values.shape[-1]
        mu = self.mu
        return -0.5 * k * tt.log(tt.sum((values - mu)**2))


class DeterminsticCalibrationFactor(pm.Continuous):
    """
    The posterior distribution for the continuous deterministic approximation
    for estimation of the fluorescence calibration factor, alpha.

    .. math::

    g(\alpha\, \vert\, k, p, [I_1, I_2]) = \alpha^{-k} \prod\limits_i^k {\Gamma \left({I_{1,i} + I_{2, i} \over \alpha} + 1\right) \over \Gamma \left({I_{1,i} \over \alpha} + 1\right)\Gamma \left({I_{2,i} \over \alpha} + 1 \right)} p^{I_{1,i} / \alpha} (1 - p)^{I_{2, i} / \alpha}


    Parameters
    ----------
    """

    def __init__(self, I_1=None, I_2=None, p=0.5, *args, **kwargs):
        super(DeterminsticCalibrationFactor, self).__init__(*args, **kwargs)
        self.p = p = pm.theanof.floatX(tt.as_tensor_variable(p))
        self.I_1 = I_1
        self.I_2 = I_2

    def logp(self, value):
        I_1 = self.I_1
        I_2 = self.I_2
        p = self.p
        k = len(I_1)

        # Compute the copy numbers with the given alpha.
        n1 = I_1 / value
        n2 = I_2 / value
        ntot = n1 + n2

        # Compute the pieces of the posterior
        binom = tt.sum(tt.gammaln(ntot + 1)) - \
            tt.sum(tt.gammaln(n1 + 1)) - tt.sum(tt.gammaln(n2 + 1))
        prob = tt.sum(n1) * tt.log(p) + tt.sum(n2) * tt.log(1 - p)

        return -k * tt.log(value) + binom + prob


class Jeffreys(pm.Continuous):
    """
    Jeffreys prior for a scale parameter.

    Parameters
    ----------
    lower : float, > 0
        Minimum value the variable can take.
    upper : float, > `lower`
        Maximum value the variable can take.
    Returns
    -------
    output : pymc3 distribution
        Distribution for Jeffreys prior.

    Notes
    -----
    This class was adopted from Justin Bois
    github.com/justinbois/bebi103
    """

    def __init__(self, lower=None, upper=None, transform='interval',
                 *args, **kwargs):
        # Check inputs
        if lower is None or upper is None:
            raise RuntimeError('`lower` and `upper` must be provided.')

        if transform == 'interval':
            transform = pm.distributions.transforms.interval(lower, upper)
        super(Jeffreys, self).__init__(transform=transform, *args, **kwargs)
        self.lower = lower = pm.theanof.floatX(tt.as_tensor_variable(lower))
        self.upper = upper = pm.theanof.floatX(tt.as_tensor_variable(upper))

        self.mean = (upper - lower) / tt.log(upper / lower)
        self.median = tt.sqrt(lower * upper)
        self.mode = lower

    def logp(self, value):
        lower = self.lower
        upper = self.upper
        return pm.distributions.dist_math.bound(
            -tt.log(tt.log(upper / lower)) - tt.log(value),
            value >= lower, value <= upper)


def ReparameterizedNormal(name=None, mu=None, sd=None, shape=1):
    """
    A reparameterized (non-centered) normal distribution. This allows for
    more efficient sampling using PyMC3

    Parameters
    ----------
    name :  string
        The name of the RV. The reparameterized version will have this name prepended with "offset_"
    mu : float
        Mean of the normal distribution.
    sd: float
        The standard deviation if the distribtion.
    shape : int
        The shape of the RV. Default is 1

    Returns
    -------
    var : PyMC3 RV object
        The reparameterized distribution.
    """
    if name is None:
        raise RuntimeError("`name` must be provided.")
    if mu is None:
        raise RuntimeError("`mu` must be provided.")
    if sd is None:
        raise RuntimeError("`sd` must be provided.")
    if type(name) is not str:
        raise TypeError(
            "expected type(name) to be string, got {0}.".format(type(name)))

    # Compute the offset.
    offset_var = pm.Normal('offset_{0}'.format(name), mu=0, sd=1, shape=shape)

    # Define the reparameterized variable.
    var = pm.Deterministic(name, mu + offset_var * sd)
    return var


def log_prior():
    """
    Computes the log prior for alpha and I_2 as needed in the model. These are taken to be uniform on the range of the camera bitdepth, 0 - 65535.
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
        The best-fit value for alpha
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
