# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np 
import pandas as pd 
import mwc.bayes

# Load the compiled data
fluct_data = pd.read_csv('../../data/compiled_fluctuations.csv')
fc_data = pd.read_csv('../../data/compiled_foldchange.csv')

# Perform appropriate background subtraction
df = []
for g, d in fluct_data.groupby(['carbon', 'date', 'run_no']):
    _fc_data = fc_data.loc[(fc_data['carbon']==g[0]) & (fc_data['date']==g[1]) & (fc_data['run_number']==g[2])]