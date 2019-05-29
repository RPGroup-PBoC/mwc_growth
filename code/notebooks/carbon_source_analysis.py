#%% [markdown]
# # Complete Analysis of Carbon Source Variation
# 
# © 2019 Griffin Chure. This work is licensed under a [Creative Commons
# Attribution License CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).
# All code contained herein is licensed under an [MIT
# license](https://opensource.org/licenses/MIT).
# 
# --- 

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.bayes
import mwc.stats
import mwc.model
import mwc.viz
import bokeh.io
import bokeh.plotting
bokeh.io.output_notebook()
colors, color_list = mwc.viz.bokeh_theme()
constants = mwc.model.load_constants()

#%% [markdown]
# In this notebook, we complete analysis of the LacI repressor titration under
# varying carbon sources. We begin by inferring the calibration factor given
# all of the data, which we model in a hierarchical fashion. With a repressor
# count in place, we then estimate the DNA binding energy of each carbon source
# and compare to that reported in Garcia & Phillips, 2011.
#%% [markdown]
# To begin, we will load the lineage information to compute the calibration factor
#%% [markdown]
# ## Computing a fluorescence calibration factor 
#%% [markdown]
# We begin by loading the data, computing the mean autofluoresnced, and defining the total integrated fluorescence for each cell 

#%%
# Load the lineage and fold-change data
lineages = pd.read_csv('../../data/compiled_fluctuations.csv')
fc_data = pd.read_csv('../../data/compiled_fold_change.csv')

# Isolate the autofluorescence data. 
auto_data = fc_data[fc_data['strain']=='auto'].copy()


# Go through and subtract teh autofluorescence information.
for g, d in lineages.groupby(['carbon', 'date', 'run_no']):
    # Extract the mean autofluorescence value for each replicate
    _auto = auto_data[(auto_data['carbon']==g[0]) &(auto_data['date']==g[1]) &
                      (auto_data['run_number']==g[2])].copy()
    lineages.loc[(lineages['carbon'] == g[0]) & (lineages['date']==g[1]) & 
                  (lineages['run_no']==g[2]), 
                  'mean_auto'] = _auto['mean_mCherry'].values.mean()

# Try subtracting the total mean autofluorescnece. 
auto = auto_data['mean_mCherry'].median()
lineages['mean_auto'] = auto

# Subtract the autofluorescence and compute the total intensity.
lineages['I1_tot'] = lineages['area_1'] * (lineages['I_1'] - lineages['mean_auto'])
lineages['I2_tot'] = lineages['area_2'] * (lineages['I_2'] - lineages['mean_auto'])
lineages['summed'] = lineages['I1_tot'] + lineages['I2_tot']
lineages['sq_diff'] = (lineages['I1_tot'] - lineages['I2_tot'])**2
lineages.dropna(inplace=True)
lineages = lineages[(lineages['I1_tot'] >= 0) & (lineages['I2_tot'] >= 0)].copy()




#%% [markdown]
# With this in hand, we can now iterate through each unique day, run_number,
# and carbon source to infer the calibration factor.
model = mwc.bayes.StanModel('../stan/hierarchical_calibration_factor.stan',
                           force_compile=True) 
stat_dfs = []
for g, d in lineages.groupby(['carbon']):
    d = d.copy()
    d['id'] = d.groupby(['date', 'run_no']).ngroup() + 1
    data_dict = {'J_exp':d['id'].max(), 'N_fluct':len(d), 'index_1':d['id'], 
                'I_1':d['I1_tot'], 'I_2':d['I2_tot']}
    samples = model.sample(data_dict)
    stats = model.summarize_parameters()
    stats['carbon'] = g
    stat_dfs.append(stats)
stat_dfs = pd.concat(stat_dfs)
stat_dfs.to_csv('../../data/calibration_factor_summary.csv', index=False)
#%%[markdown]
# Now, let's build a Bokeh app to look at each day and it's respective
# calibration factor fitting. 
#%%
# Instantiate the figure. 
# carbon = 'glycerol'
# carbon = 'glucose'
carbon = 'acetate'
ax = bokeh.plotting.figure(width=600, height=400,
                           x_axis_type='log',
                           y_axis_type='log',
                           x_axis_label='sum total fluorescence [a.u.]',
                           y_axis_label='squared fluctuations')

c_data = lineages[lineages['carbon']==carbon]
I_tot = np.logspace(1, 6, 100)

iter = 0 
for g, d in c_data.groupby(['date', 'run_no']):
    ax.circle(d['summed'], d['sq_diff'], size=1, color=color_list[iter])
    iter +=1

#     low = param['hpd_min'].values[0] * I_tot
#     high = param['hpd_max'].values[0] * I_tot
#     mwc.viz.fill_between(ax, I_tot, low, high, color=color_list[iter], alpha=0.5)
    # iter += 1

stats = stat_dfs[(stat_dfs['carbon']==carbon) &\
     (stat_dfs['parameter']=='alpha_1')][['hpd_min', 'hpd_max']].values[0]
low = stats[0] * I_tot
high = stats[1] * I_tot#
mwc.viz.fill_between(ax, I_tot, low, high, color='dodgerblue', alpha=0.5)
ax.legend.location = 'top_left'
bokeh.io.show(ax)

# %%[markdown]
# That looks okay, it's always hard to tell what is a good "fit" with this model
# as the variance of the binomial is often so large. The credible regions look
# believable and are about the same for every carbon source. Using the credible
# reigons, we can compute the lower and upper bound for the repressor copy
# number given the total fluroescence. Note that here I will need to subtract
# the autofluorescence.
#%%
# Subtract the autofluorescence from each sample and compute the total
# intensity.  This is loppy but it works for now.
valid_dfs = []
for g, d in fc_data.groupby(['carbon', 'date', 'run_number']):
    _auto = d[d['strain']=='auto']
    _delta = d[d['strain']=='delta']

    # Isolate the calibration factor credible regions. 
    _stats = stat_dfs[(stat_dfs['carbon']==g[0]) & (stat_dfs['parameter']=='alpha_1')]
    if len(_stats) > 0:
        d = d[['carbon', 'date', 'run_number', 'mean_mCherry', 'mean_yfp', 'area_pix',
               'atc_ngml', 'strain']].copy()
        d
        d['auto_mCherry'] = _auto['mean_mCherry'].mean()
        d['auto_yfp'] = _auto['mean_yfp'].mean()
        d['delta_yfp'] = np.mean(_delta['area_pix'].values * (_delta['mean_yfp'] - _auto['mean_yfp'].mean()))
        d['Itot_mCherry'] = d['area_pix'].values * (d['mean_mCherry'].values - d['auto_mCherry'].values)
        d['Itot_yfp'] = d['area_pix'].values * (d['mean_yfp'].values - d['auto_yfp'])
        d['alpha_min'] = _stats['hpd_min'].values[0]
        d['alpha_max'] = _stats['hpd_max'].values[0]
        d['alpha_median'] = _stats['median'].values[0]
        d['rep_min'] = d['Itot_mCherry'].values / _stats['hpd_max'].values[0]
        d['rep_max'] = d['Itot_mCherry'].values / _stats['hpd_min'].values[0]
        d['rep_median'] = d['Itot_mCherry'].values /  _stats['median'].values[0]
        d['fold_change'] = d['Itot_yfp'].values / d['delta_yfp'].values
        valid_dfs.append(d)
valid_dfs = pd.concat(valid_dfs)
valid_dfs.to_csv('../../data/folchange_represors.csv', index='False')
#%%
rep_range = np.logspace(0, 3, 200)
theo = (1 + 0.99 * (rep_range / 4.6E6) * np.exp(13.9))**-1
# Plot the clouds of fold-change.
ax = bokeh.plotting.figure(width=500, height=300, x_axis_type='log',
                           y_axis_type='log',
                           x_axis_label='repressors per cell',
                           y_axis_label='fold-change',
                           x_range=[1, 1000])

# Plot the clouds. 
c_color = {'glucose':colors['purple'], 
           'acetate':colors['orange'], 
           'glycerol':colors['blue']}

for g, d in valid_dfs.groupby(['carbon']):
    if g == 'acetate':
        d = d[d['strain']=='dilution']
        ax.circle(d['rep_median'], d['fold_change'], color=c_color[g], size=0.75,
                legend=g, alpha=0.5)

ax.line(rep_range, theo, color='black', legend='theory')

ax.legend.location="bottom_left"
bokeh.io.show(ax)

valid_dfs['alpha_min'].unique()
#%% [markdown]
# Interesting. It looks like *all* samples are more repressed than expected,
# including glucose. Playing around with the binding energy makes me think that
# this O2 sequence is actually O1. I will do some sequencing to figure it all
# out.  In the mean time, let's bin this data a bit to get a better sense of how
# the means track the theory. 

#%% [markdown]
# ## Inferring the DNA binding energy



#%%
model = mwc.bayes.StanModel('../stan/single_cell_epRA.stan', force_compile=True)

samples_dfs = []
stat_dfs = []

for g, d in valid_dfs.groupby(['carbon']):
    print(f'Sampling {g}...')
    data_dict = {'N':len(d), 'R':d['rep_median'], 'fc':d['fold_change']}
    samples = model.sample(data_dict, iter=5000, control=dict(adapt_delta=0.99))
    _df = samples[1]
    _df['carbon'] = g 
    _stats = model.summarize_parameters()
    _stats['carbon'] = g
    stat_dfs.append(_stats)
    samples_dfs.append(_df)
    print('... finished!')
#%%
stats = pd.concat(stat_dfs)


#%%
