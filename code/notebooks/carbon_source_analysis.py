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
fc_data.dropna(inplace=True)

# Isolate the autofluorescence data. 
auto_data = fc_data[fc_data['strain']=='auto'].copy()

# Go through and subtract teh autofluorescence information.
for g, d in lineages.groupby(['carbon', 'date', 'run_no']):
    # Extract the mean autofluorescence value for each replicate
    _auto = auto_data[(auto_data['carbon']==g[0]) &(auto_data['date']==g[1]) &
                      (auto_data['run_number']==g[2])].copy()
    fc_data.loc[(fc_data['carbon']==g[0]) & (fc_data['date']==g[1]) & (fc_data['run_number']==g[2]), 'mean_auto'] = (_auto['mean_mCherry'].values).mean()
    lineages.loc[(lineages['carbon'] == g[0]) & (lineages['date']==g[1]) & 
                  (lineages['run_no']==g[2]), 
                  'mean_auto'] = (_auto['mean_mCherry'].values).mean()

for g, d in lineages.groupby(['carbon']):
        # Try subtracting the total mean autofluorescnece. 
        auto = d['mean_auto'].mean()
        lineages.loc[lineages['carbon']==g, 'mean_auto_c'] = auto


# Subtract the autofluorescence and compute the total intensity.
lineages['I1_tot'] = lineages['area_1'] * (lineages['I_1'] - lineages['mean_auto_c'])
lineages['I2_tot'] = lineages['area_2'] * (lineages['I_2'] - lineages['mean_auto_c'])
lineages['summed'] = lineages['I1_tot'] + lineages['I2_tot']
lineages['fluct'] = (lineages['I1_tot'] - lineages['I2_tot'])**2
lineages.dropna(inplace=True)
lineages = lineages[(lineages['I1_tot'] >= 0) & (lineages['I2_tot'] >= 0)].copy()


#%% [markdown]
# With this in hand, we can now iterate through each unique day, run_number,
# and carbon source to infer the calibration factor.
# model = mwc.bayes.StanModel('../stan/calibration_factor.stan')

df = pd.DataFrame([])
for g, d in lineages.groupby(['carbon']):
    d = d.copy()
    opt, std = mwc.bayes.estimate_calibration_factor(d['I1_tot'], d['I2_tot'])
    df = df.append({'carbon':g, 'alpha_opt':opt, 'alpha_std':std}, ignore_index=True)
   
#     data_dict = {'J_exp':d['id'].max(), 'N':len(d), 'index_1':d['id'], 
                # 'I1':d['I1_tot'], 'I2':d['I2_tot']}
#     samples = model.sample(data_dict)
#     stats = model.summarize_parameters()
#     stats['carbon'] = g
#     stat_dfs.append(stats)

# stat_dfs = pd.concat(stat_dfs)
# stat_dfs.to_csv('../../data/calibration_factor_summary.csv', index=False)

# %% [markdown]
# Note that here we took a pooled approach in estimation of the calibration factor
#%%
# Instantiate the figure. 
carbon = 'acetate'
# carbon = 'glucose'
# carbon = 'glycerol'
ax = bokeh.plotting.figure(width=600, height=400,
                           x_axis_type='log',
                           y_axis_type='log',
                           x_axis_label='sum total fluorescence [a.u.]',
                           y_axis_label='squared fluctuations')

c_data = lineages[lineages['carbon']==carbon]
I_tot = np.logspace(1, 6, 100)

iter = 0 
for g, d in c_data.groupby(['date', 'run_no']):
    # Compute bins
    ax.circle(d['summed'], d['fluct'], size=1, color='black')
    iter +=1


bins = mwc.stats.bin_by_events(c_data, bin_size=200)
ax.circle(bins['summed'], bins['fluct'], color='tomato', size=5) 


# stats = stat_dfs[(stat_dfs['carbon']==carbon) &\
#      (stat_dfs['parameter']=='alpha')][['hpd_min', 'hpd_max']].values[0]

opt = df[df['carbon']==carbon]['alpha_opt'].values[0] * I_tot
low = stats[0] * I_tot
high = stats[1] * I_tot
ax.line(I_tot, opt, color='tomato', alpha=0.75, line_width=4)
# mwc.viz.fill_between(ax, I_tot, low, high, color='tomato', alpha=0.5)

ax.legend.location = 'top_left'
bokeh.io.show(ax)

# %%[markdown]
# ## Exploring the limits of measurement 
# Before we can calculate the repressors vs fold-change, we need to konw what
# our noise floor is -- i.e. how many repressors per cell are calculated for the
# autofluorescence samples? Let's look at a distribution to figure it out. 

#%%
_auto_data = fc_data[fc_data['strain']=='auto']
for g, d in _auto_data.groupby('carbon'):
        # _stats = stat_dfs[(stat_dfs['carbon']==g) &\
                        # (stat_dfs['parameter']=='alpha')]['median']
        opt = df[df['carbon']==g]['alpha_opt'].values[0]
        _mean_auto = (d['mean_mCherry']).mean()
        I_tot = d['area_pix'] * (d['mean_mCherry'] - _mean_auto)
        rep_mean = I_tot / opt
        _auto_data.loc[_auto_data['carbon']==g, 'mean_rep'] = rep_mean

#%%
ax = bokeh.plotting.figure()
_colors = {'glucose':colors['purple'], 'glycerol':colors['blue'], 
                'acetate':colors['orange']}
for g, d in _auto_data.groupby('carbon'):
        x, y = np.sort(d['mean_rep']), np.linspace(0, 1, len(d))
        ax.step(x, y, color=_colors[g], line_width=1)
        ax.circle(x, y, color=_colors[g], legend = g)

ax.legend.click_policy = 'hide'
bokeh.io.show(ax)

#%% [markup]
# Looking at the ECDFs, it seems like ~ 20 repressors per cell is the maximum
# observed count for hte autofluorescence. I will use this as my lower bound and
# drop any calculated repressors below this floor.  

#%%
lower_bound = 5 

# Compute the repressor count, only go with the median. 
# dilution_data = fc_data[fc_data['strain']=='dilution']
_fc_data = fc_data.copy()
for g, d in _fc_data.groupby(['carbon', 'strain']):
        # Isolate the calibration factor median
        # _stats = stat_dfs[(stat_dfs['carbon']==g[0]) & 
                        #   (stat_dfs['parameter']=='alpha')]['median']
        opt  = df[df['carbon']==g[0]]['alpha_opt'].values[0]
        _fc_data.loc[(_fc_data['carbon']==g[0]) & 
                        (_fc_data['strain']==g[1]), 'repressors'] = np.round(d['area_pix'] * (d['mean_mCherry'] - d['mean_auto']) / opt)
        
_fc_data['repressors'] = _fc_data['repressors'].values.astype(int)
_fc_data = _fc_data[_fc_data['repressors'] > lower_bound]

# %% [markdown]
# Now, let's separate each carbon source into some bins, compute the means, and
# plot it. 
#%%
ax = bokeh.plotting.figure(x_axis_type='log', y_axis_type='log')
_colors = {'glucose':colors['purple'], 'glycerol':colors['light_grey'], 'acetate':colors['orange']}
for g, d in _fc_data.groupby(['carbon']):
        delta  = d[d['strain']=='delta'].copy()
        dil = d[d['strain']=='dilution'].copy()
        mean_auto = d[d['strain']=='auto']['mean_yfp'].mean()
        mean_delta = delta['mean_yfp'] - mean_auto
        dil['fc'] = (d['mean_yfp'] - mean_auto) / mean_delta
        
        # Separate into bins. 
        dil.sort_values('repressors', inplace=True)
        # _bins = mwc.stats.bin_by_events(dil, 400, sortby='repressors', average=['fold_change', 'repressors'])
        # ax.circle(_bins['repressors'], _bins['fold_change'], color=_colors[g], legend=g)
        ax.circle(dil['repressors'], dil['fold_change'], color=_colors[g], legend=g,
        size=2, alpha=0.4)

# Plot the theory
rep_range = np.logspace(0, 3.5)
ax.line(rep_range, (1 + (rep_range / 4.6E6) * np.exp(13.9))**-1, color='black', legend='theory')
ax.legend.click_policy = 'hide'
bokeh.io.show(ax)
#%%[markdown]
# ## Inference of Fold-Change
# I'm going balls-to-the-wall for inference, so let's just do it for the
# fold-change as well!


#%%
model = mwc.bayes.StanModel('../stan/hierarchical_foldchange.stan', force_compile=True)




#%%
#%%
import tqdm
fc_dfs = []
for g, d in _fc_data.groupby(['carbon']):
        # Separate by strain
        auto_df = d[d['strain']=='auto']
        delta_df = d[d['strain']=='delta']
        dil_df = d[d['strain']=='dilution']

        # Drop repressors below the threshold
        dil_df = dil_df[dil_df['repressors'] > lower_bound]

        dil_df['bin'] = pd.cut(dil_df['repressors'], 25)

        # Group by the bins 
        _stats_dfs = []

        for _g, _d in tqdm.tqdm(dil_df.groupby(['bin']), desc='Inferring fold-change'):
                # Assemble the data dictionary
                data_dict = {'J':1, 'idx':np.ones(len(_d)).astype(int),
                     'N_dilution':len(_d), 'N_auto':len(auto_df), 
                     'N_delta':len(delta_df), 'auto_yfp':auto_df['mean_yfp'],
                     'delta_yfp':delta_df['mean_yfp'],
                     'dilution_yfp':_d['mean_yfp']}

                # Sample
                fit, samples = model.sample(data_dict)
                stats = model.summarize_parameters()
                stats['mean_rep'] = np.mean(_d['repressors'])
                stats['carbon']=g
                stats['repressors'] = _g
                fc_dfs.append(stats)

#%%
