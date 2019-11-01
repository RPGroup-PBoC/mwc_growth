#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
colors, color_list = mwc.viz.personal_style()

# %%
plates = pd.read_csv('../../data/compiled_growth_plates.csv')
stats = pd.read_csv('../../data/compiled_growth_statistics.csv')
plates['temp_C'] = np.round(plates['temp_C'])
plates['time_min'] = np.round(plates['time_min'])
# plates = plates[plates['time_min'] <= 600]

# Load specific examples of GP processing. 
gp_output = pd.read_csv('../processing/growth_curves/20190206_r1_37C_mixedmedia_O2_growth/output/delta_glycerol/gp_output.csv')
plate = pd.read_csv('../processing/growth_curves/20190206_r1_37C_mixedmedia_O2_growth/output/growth_plate.csv')
# %%
plate

# %%
fig, ax = plt.subplots(1, 1, dpi=100)
ax.plot(gp_output['time'], gp_output['OD_raw_data'], 'k,')

# %%
