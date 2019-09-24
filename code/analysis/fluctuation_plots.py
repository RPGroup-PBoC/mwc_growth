#%%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import bokeh.layouts
from bokeh.transform import factor_cmap
import bokeh.palettes
import mwc.bayes
import mwc.viz
import mwc.stats
bokeh.io.output_notebook()
colors, color_list = mwc.viz.bokeh_theme()

# %% Load the fluctuation data 
fluct_data = pd.read_csv('../../data/analyzed_fluctuations.csv')
foldchange = pd.read_csv('../../data/analyzed_foldchange.csv')

#%%
CARBON = 'glucose'
TEMP = 37

#Isolate the fluctuations
flucts = fluct_data[(fluct_data['carbon']==CARBON) & 
                    (fluct_data['temp']==TEMP)]
fc = foldchange[(foldchange['carbon']==CARBON) & 
                (foldchange['temp']==TEMP) &
                (foldchange['strain']=='dilution')]

summarized = fc.groupby(['atc_ngml']).mean().reset_index()
# Assign replicate IDs
flucts['idx'] = flucts.groupby(['date', 'run_no']).ngroup() + 1

# Set up the figure
fluct_fig = bokeh.plotting.figure(width=600, height=400, x_axis_type='log',
                            y_axis_type='log', 
                            x_axis_label='summed intensities [a.u.]',
                            y_axis_label='square fluctuations [a.u.]')
fc_fig = bokeh.plotting.figure(width=600, height=400, x_axis_type='log',
                            y_axis_type='log', 
                            x_axis_label='repressors per cell',
                            y_axis_label='fold-change')
R = np.logspace(0, 3, 500)
fc_fig.line(R, (1 + (R/5E6) * np.exp(13.9))**-1, color='tomato')
fc_fig.circle('repressors', 'fold_change', source=summarized, color='black', size=1)

fluct_fig.circle('summed', 'fluct', source=flucts, size=2)
lay = bokeh.layouts.column(fc_fig, fluct_fig)
bokeh.io.show(lay)


#%%
