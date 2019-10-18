#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
colors, color_list = mwc.viz.personal_style()


# Define some experimental constants. 
epRA_ref = -13.9 # in kBT
epRA_sig = 0.2 # in kBT from Garcia & Phillips 2011
epAI_ref = 4.5 # in kBT from Razo-Mejia et al. 2018
temp_ref = 37
temp_range = np.linspace(30, 45, 100) + 273.15  # in kBT

# Load the statistics for the inferred binding energies. 
fit_params = pd.read_csv('../../data/DNA_binding_energy_summary.csv')
fit_params = fit_params[fit_params['carbon']=='glucose']

#%%
# Instantiate the figure
fig, ax = plt.subplots(1, 1, figsize=(3.42, 3.42), dpi=150)


# Plot the theoretical predictions
simple_min = (epRA_ref - epRA_sig) * (temp_ref + 273.15) / (temp_range)
simple_max = (epRA_ref + epRA_sig) * (temp_ref + 273.15) / (temp_range)

ent_min = 
ax.fill_between(temp_range, simple_min, simple_max, color='k', alpha=0.25)

# Plot the inferred binding energies
for g, d in fit_params.groupby('temp'):
    RA_median, RA_hpd_min, RA_hpd_max = d[d['parameter']=='epRA'][
                                ['median', 'hpd_min', 'hpd_max']].values[0]
    ax.plot(g + 273.15, RA_median, 's')


#%%



#%%
