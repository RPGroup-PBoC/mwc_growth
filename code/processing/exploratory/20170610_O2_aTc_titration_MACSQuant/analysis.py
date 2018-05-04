import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import glob
import sys
sys.path.insert(0, '../../../')
import mwc_growth as mwc
mwc.set_plotting_style()

# Define the constants.
fname = glob.glob('output/*.csv')

# Load the data.
data = pd.read_csv(fname[0])
lac_data = data[data['rbs']=='lacI']
plt.plot(lac_data['ATC'], lac_data['mean_mCherry_H']/lac_data['mean_mCherry_H'].loc[8], 'r--o', label='mCherry')
plt.plot(lac_data['ATC'], lac_data['mean_YFP_H'] / lac_data['mean_YFP_H'].loc[2], 'b--o', label='YFP')
plt.xlabel('anhydro-tetracycline (ng / mL)')
plt.ylabel('fractional mean fluorescence of max')
plt.legend()
plt.savefig('output/20170610_initial_aTc_titration.png', bbox_inches='tight')
