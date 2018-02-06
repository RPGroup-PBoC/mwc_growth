# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import statsmodels.tools.numdiff as smnd
import scipy.optimize
import glob


def ecdf(data):
    """
    Computes the empirical cumulative distribution function for a collection
    of provided data.
    """
    return np.sort(data), np.arange(0, len(data)) / len(data)
