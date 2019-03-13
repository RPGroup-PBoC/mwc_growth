import numpy as np
import os
import pandas as pd

filename = '20190102_r2_37C_glucose_O2_foldchange.csv'
_fc = pd.read_csv(filename)
fc_pruned = _fc[_fc['atc_ngml']!=0.4]
fc_pruned.to_csv(filename, index=False)