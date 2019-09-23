#%%
import numpy as np
import bokeh.io
import bokeh.plotting
import mwc.bayes
import mwc.viz
import mwc.stats

# %% Load the fluctuation data 
fluct_data = pd.read_csv