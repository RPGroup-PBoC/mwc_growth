# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import pickle
import pystan
import pandas as pd
import scipy.special
import scipy.optimize
import statsmodels.tools.numdiff as smnd
from .viz import bokeh_traceplot


class StanModel(object):
    R"""
    Custom StanModel class for crafting and sampling from Stan
    models.
    """
    def __init__(self, model, data_dict, samples=None, force_compile=False):
        """
        Parameters
        ----------
        model: str
            Relative path to saved Stan model code. To deter bad habits,
            this class does not accept a string as the model code. Save 
            your Stan models. 
        data_dict: dictionary
            Dictonary of all data block parameters for the model.
        force_compile: bool
            If True, model will be forced to compile. If False, 
            a precompiled file will be loaded if present. 
        """
        self.model = loadStanModel(model, force=force_compile)
        self.data = data_dict
        self.samples = samples
        self.df = None
        
    def sample(self, iter=2000, chains=4, return_df=True, **kwargs):
        """
        Samples the assembled model given the supplied data dictionary
        and returns output as a dataframe.
        """
        self.chains = chains
        self.iter = iter
        print('Beginning sampling...')
        self.samples = self.model.sampling(self.data, 
                        chains=chains, iter=iter, **kwargs)
        print('finished sampling!')
        if return_df:
            self.df = self.samples.to_dataframe(diagnostics=True)
            return [self.samples, self.df]
        else:
            return self.samples
    
    def check_divergence(self, thresh=0.0005, return_values=True, quiet=False):
        """
        Computes the fraction of diverging samples
        """
        if self.samples == None:
            raise RuntimeError('Divergence is not defined without sampling. Please sample your model first')
        if self.df == None:
            self.df = self.samples.to_datframe(diagnostics=True)

        n_div = np.sum(self.df['divergent__']) 
        div_frac = n_div / len(self.df) * 100
        if div_frac == 0:
            statement = "No diverging samples found. Nicely done."
        if div_frac < thresh:
            statement = "Diverging samples below {} % ({} of {} samples diverging).".format(thresh * 100, n_div, len(self.df)) 
        else:
            statement = "Warning, {} % of samples are diverging. Reparameterize your model or adjust adapt_delta above 0.8.".format(div_frac)
        
        if quiet is not False:
            print(statement)
        if return_values:
            return {'statement':statement, 'n_diverging':n_div, 'n_samples': len(self.df), 'diverging_fraction':div_frac}
        else:
            return statement

    def check_rhat(self, return_values=True, quiet=False):
        """
        Determines the Gelman-Rubin statistic (R-hat). If 0.9 < r-hat < 1.1, the sampler has converged. 
        """
        if self.samples == None:
            raise RuntimeError('R-hat not defined without sampling. Please sample your model first.')
        if self.df == None:
            self.df = self.samples.to_dataframe(diagnostics=True)
        raise RuntimeError('Not yet implemented!')

    def check_n_effective(self, thresh=0.001, return_values=True, quiet=False):
        if self.samples == None: 
           raise RuntimeError('n_effective / N not defined without sampling. Please sample your model first.')
        if self.df == None:
            self.df = self.samples.to_dataframe(diagnostics=True)
        raise RuntimeError('Not yet implemented!')

    def check_diagnostics(self, save_summary=False, fname=None, return_values=True, quiet=False):
        """
        Checks all sampling diagnostics. 

        Parameters
        ----------
        save_summary: bool
            If True, a summary file will be saved. fname is required. 
        fname: str
            Desired filename of summary file. Only required if save_summary is True.
        return_values: bool
            If True, a dictionary of diagnostics is returned.
        quiet: bool
            If True, summary will not be printed to screen. Default is False.
        """
        raise RuntimeError('Not yet implemented!')
        
    def traceplot(self, varnames=None):
        """
        Shows the sampling trace and distributions for desired varnames
        See documentation for mwc.viz.bokeh_traceplot for more details.
        """
        return bokeh_traceplot(self.samples, varnames=varnames)

    

    
    

def loadStanModel(fname, force=False):
    """Loads a precompiled Stan model. If no compiled model is found, one will be saved."""
    # Identify the model name and directory structure
    rel, sm_dir = fname.split('/stan/')
    sm_name = sm_dir.split('.stan')[0]
    pkl_name = f'{rel}/stan/{sm_name}.pkl' 
    # Check if the model is precompiled
    if (os.path.exists(pkl_name)==True) and (force != True):
        print('Found precompiled model. Loading...')
        model = pickle.load(open(pkl_name, 'rb'))
        print('finished!')
    else:
        print('Precompiled model not found. Compiling model...')
        model = pystan.StanModel(fname)
        print('finished!')
        with open(pkl_name, 'wb') as f:
            pickle.dump(model, f)      
    return model
    

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

    # Assemble the log posterior.
    logpost = change_of_var + binom + prob
    return prefactor * (logpost)


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

   