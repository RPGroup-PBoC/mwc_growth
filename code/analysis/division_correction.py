# -*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import bokeh.palettes
import mwc.viz
import mwc.bayes
import mwc.stats
colors, color_list = mwc.viz.bokeh_theme()
bokeh.io.output_notebook()
# %%
lineages = pd.read_csv('../../data/compiled_fluctuations.csv')
data = pd.read_csv('../../data/compiled_fold_change.csv')

# %% Look at the distributions of areas
x = np.sort(data['area_pix'] * 0.065**2)
y = np.arange(0, len(x), 1) / len(x)

ax = bokeh.plotting.figure(width=600, height=400,
                            x_axis_label='area [square µm]',
                            y_axis_label='cdf')

iter = 0
for g, d in data[data['carbon']=='glucose'].groupby(['atc_ngml']):
    x = np.sort(d['area_pix'] * 0.065**2)
    y = np.arange(0, len(x), 1) / len(x)
    ax.step(x, y, color=color_list[iter], legend=str(g), line_width=2)
    iter += 1

ax.legend.click_policy = 'hide'
bokeh.io.show(ax)
#%%
ax = bokeh.plotting.figure()
iter = 0
for g, d in data.groupby(['atc_ngml', 'carbon', 'temp', 'date', 'run_number']):
    ax.circle(g[0], d['area_pix'].mean() * 0.065**2)
    # ax.line(g, d['area_pix'].mean() * 0.065**2)

bokeh.io.show(ax)

#%%
