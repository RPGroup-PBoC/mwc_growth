# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import mwc.io
import mwc.stats
import mwc.bayes

# Load the various data. 
fc_data = pd.read_csv('../../data/compiled_fold_change.csv')
growth_data = pd.read_csv('../../data/mean_area_growth_samples.csv')

