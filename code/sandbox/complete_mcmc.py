# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pystan 
sys.path.insert(0, '../../')
import mwc.viz
import mwc.stats
import mwc.bayes
import scipy.special
colors = mwc.viz.personal_style()

# Generate a fake dataset.
alpha = 150
n_tot = np.random.randint(1, 1E3, 200)
n1 = np.random.binomial(n_tot, 0.5)
n2 = n_tot - n1
I1 = alpha * n1 + np.random.normal(0, 10, len(n1))
I2 = alpha * n2 + np.random.normal(0, 10, len(n1))

# %%
# Load the stan model
model = pystan.StanModel('complete_mcmc.stan')

#%% Set up the data dictionary.
data_dict = {'N': len(n1), 'I1':I1, 'I2':I2, 'p':0.5}
samples = model.sampling(data_dict, iter=10000, chains=2, thin=10)
samples