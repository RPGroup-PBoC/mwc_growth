# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import bebi103.viz
import mwc.bayes
import mwc.viz

# Load in the compiled data
fluct_data = pd.read_csv('../../data/compiled_fluctuations.csv')
fc_data = pd.read_csv('../../data/compiled_fold_change.csv')

# Isolate a single carbon source for dnow
fluct_data = fluct_data[fluct_data['carbon']=='glucose']
fc_data = fc_data[fc_data['carbon']=='glucose']


# Determine the mean autofluorescence. 
auto = fc_data[fc_data['strain']=='auto']
mean_auto = np.mean(auto['mean_mCherry'] - auto['mCherry_bg_val'])

# Perform background subtraction for fluctuation data. 
fluct_data['I_1'] = (fluct_data['I_1'] - (fluct_data['bg_val'] + mean_auto)) * fluct_data['area_1']
fluct_data['I_2'] = (fluct_data['I_2'] - (fluct_data['bg_val'] + mean_auto)) * fluct_data['area_2']
fc_data['I_tot'] = (fc_data['mean_mCherry'] - (fc_data['mCherry_bg_val'] + mean_auto)) * fc_data['area_pix']

# Assign various identifiers
fluct_data['run_idx'] = fluct_data.groupby(['date', 'run_no']).ngroup() + 1
dil_data = fc_data[fc_data['strain']=='dilution']
dil_data['conc_idx'] = dil_data.groupby(['atc_ngml']).ngroup() + 1
dil_data['run_idx'] = dil_data.groupby(['atc_ngml', 'date', 'run_number']).ngroup() + 1


# Generate the data dictionary. 
data_dict = {'J_exp': np.max(fluct_data['run_idx']),
             'N_fluct': len(fluct_data),
             'index_1': fluct_data['run_idx'],
             'J_conc': np.max(dil_data['conc_idx']),
             'J_conc_exp': np.max(dil_data['run_idx']),
             'N_mch': len(dil_data),
             'index_2': dil_data['conc_idx'],
             'index_3': dil_data['run_idx'],
             'I_1': fluct_data['I_1'],
             'I_2': fluct_data['I_2'],
             'mcherry': fc_data['I_tot']}
# Compile the model and sample. 
model = mwc.bayes.StanModel('../stan/hierarchical_calibration_factor.stan', data_dict, force_compile=True)

