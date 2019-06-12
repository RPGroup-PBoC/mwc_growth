# %%
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import mwc.viz
import mwc.io
import imp
import mwc.bayes
imp.reload(mwc.bayes)
import bokeh.io
import bokeh.plotting
import bokeh.palettes
import bokeh.transform
#%% 
bokeh.io.output_notebook()
colors, color_list = mwc.viz.bokeh_theme()

# Load the data sets
lineages = pd.read_csv('../../data/compiled_fluctuations.csv')
data = pd.read_csv('../../data/compiled_fold_change.csv')

# Keep only glucose
lineages = lineages[lineages['carbon']=='glucose'].copy()
data = data[data['carbon']=='glucose'].copy()

# Compute the mCherry autofluorescence background. 
auto_37 = data[(data['temp']==37) & (data['strain']=='delta')]['mean_mCherry'].mean()
auto_42 = data[(data['temp']==42) & (data['strain']=='delta')]['mean_mCherry'].mean()
mean_bg = lineages['bg_val'].mean()
bg_yfp = data[data['strain']!='dilution']['yfp_bg_val'].mean()
auto_dict = {37:auto_37 - mean_bg, 42:auto_42 - mean_bg}

#%%

# Load the stan model
model = mwc.bayes.StanModel('../stan/calibration_factor.stan')
dfs = []
data_dfs = []
for g, d in lineages.groupby('temp'):
    print(g)
    d = d.copy()
    auto = auto_dict[g]

    # Compute the integrated intensities
    d['I1'] = d['area_1'] * (d['I_1'] - auto - mean_bg)
    d['I2'] = d['area_2'] * (d['I_2'] - auto - mean_bg)
    d = d[(d['I1'] >= 0) & (d['I2'] >= 0)]
    data_dfs.append(d)
    # Sample the model
    _, samples = model.sample({'N':  len(d), 'I1':d['I1'], 'I2':d['I2']})
    summary = model.summarize_parameters()
    summary['temp'] = g
    dfs.append(summary)
stats = pd.concat(dfs)
data_df = pd.concat(data_dfs)


#%%
# Compute the sum and fluct for each
data_df['sum'] = data_df['I1'] + data_df['I2']
data_df['fluct'] = (data_df['I1'] - data_df['I2'])**2
d37 = data_df[data_df['temp']==37]
d42 = data_df[data_df['temp']==42]

p = bokeh.plotting.figure(x_axis_type='log', y_axis_type='log')
p.circle(x='sum', y='fluct', source=d37, color=color_list[0], legend=str(37), alpha=0.5)
p.circle(x='sum', y='fluct', source=d42, color=color_list[1], legend=str(42))
bokeh.io.show(p)

#%%


#%%
# Set up the axis to plot the fold-change data
ax = bokeh.plotting.figure(x_axis_type='log', y_axis_type='log', x_range=[1, 500])

# Compute the theory
R = np.logspace(0, 4,  200)
theory = (1 + (1 + np.exp(-4.5))**-1 * (R / 4E6) * np.exp(13.9))**-1

color_idx = {37:color_list[0], 42:color_list[1]}
for g, d in data.groupby('temp'):
    alpha=stats[(stats['temp']==g) & (stats['parameter']=='alpha')]

    # Compute the integrated mCherry and integrated YFP
    d = d.copy()
    auto = d[d['strain']=='auto']
    delta = d[d['strain']=='delta']
    auto_yfp = auto['mean_yfp'].mean() - bg_yfp
    auto_mch = delta['mean_mCherry'].mean() - mean_bg

    # COmpute the quantities.
    d['repressors'] = d['area_pix'] * (d['mean_mCherry'] - mean_bg - auto_mch) / alpha['median'].values[0]
    d['repressors_min'] = d['area_pix'] * (d['mean_mCherry'] - mean_bg - auto_mch) / alpha['hpd_min'].values[0]
    d['repressors_max'] = d['area_pix'] * (d['mean_mCherry'] - mean_bg - auto_mch) / alpha['hpd_max'].values[0]
    d['int_yfp'] = d['area_pix'] * (d['mean_yfp'] - auto_yfp - bg_yfp)
    d['fold_change'] = d['int_yfp'] / d[d['strain']=='delta']['int_yfp'].mean()

    # Plot the cloud
    d = d[(d['strain'] == 'dilution')]
    d['repressors'] = np.round(d['repressors'].values)
    d = d[(d['repressors'] >= 10) & (d['repressors'] <= 200)]

    if g == 37:
        size=2
        alpha=0.5
    else:
        size=5
        alpha=1

    _d = d.groupby(['repressors'])[['repressors_min', 'repressors_max', 'fold_change'].agg(('mean', 'sem')).reset_index()
    ax.circle(2 * _d['repressors'], _d['mean'], 
                size=size, alpha=alpha, legend=str(g), color=color_idx[g])
    ax.segment(x0=2 * _d['repressors'], x1=2 * _d['repressors'],  y0=_d['mean'] - _d['sem'], 
                y1=_d['mean'] + _d['sem'], line_width=2, 
                alpha=alpha, legend=str(g), color=color_idx[g])

ax.line(R, theory, color='black', legend='theory')
ax.legend.click_policy = 'hide'
bokeh.io.show(ax)

#%%


#%%
