# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pymc3 as pm
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


class DeterminsticCalibrationFactor(pm.PositiveContinuous):
    """
    The posterior distribution for the continuous deterministic approximation
    for estimation of the fluorescence calibration factor, alpha.

    .. math::

    g(\alpha\, \vert\, k, p, [I_1, I_2]) = \alpha^{-k} \prod\limits_i^k {\Gamma \left({I_{1,i} + I_{2, i} \over \alpha} + 1\right) \over \Gamma \left({I_{1,i} \over \alpha} + 1\right)\Gamma \left({I_{2,i} \over \alpha} + 1 \right)} p^{I_{1,i} / \alpha} (1 - p)^{I_{2, i} / \alpha}


    Parameters
    ----------
    """

    def __init__(self, alpha=None, p=0.5, *args, **kwargs):
        if alpha is None:
            raise RuntimeError('variable alpha must be supplied')
        super(DeterminsticCalibrationFactor, self).__init__(*args, **kwargs)
        self.alpha = alpha = pm.theanof.floatX(tt.as_tensor_variable(alpha))
        self.p = p = pm.theanof.floatX(tt.as_tensor_variable(p))

    def logp(self, value):
        I_1, I_2 = value
        p = self.p
        alpha = self.alpha
        k = len(I_1)

        # Compute the copy numbers with the given alpha.
        n1 = I_1 / alpha
        n2 = I_2 / alpha
        ntot = n1 + n2

        # Compute the pieces of the posterior
        binom = tt.sum(tt.gammaln(ntot + 1)) - \
            tt.sum(tt.gamaln(n1 + 1)) - tt.sum(tt.gammaln(n2 + 1))
        prob = tt.sum(n1) * tt.log(p) + tt.sum(ntot) * tt.log(1 - p)

        return -k * tt.log(alpha) + binom + prob


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


# Reparameterization Functions
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
