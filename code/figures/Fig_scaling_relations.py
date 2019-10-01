#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import mwc.stats
colors, color_list = mwc.viz.personal_style()

#%%
# Load the fluorescence data. 
data = pd.read_csv('../../data/analyzed_foldchange.csv')

# Load the growth rate data. 
stats = pd.read_csv('../../data/compiled_growth_statistics.csv')

# Properly format the growth stats
stats = stats[((stats['carbon']=='glucose') | (stats['carbon']=='acetate') |
                (stats['carbon']=='glycerol')) & 
                ((stats['temp']==37) |  (stats['temp']==32) | 
                (stats['temp']==42))] 

tidy_stats = pd.DataFrame([])
for g, d in stats.groupby(['date', 'carbon', 'temp', 'run_number']):
    growth_rate = d[d['parameter']=='max df']['value'].values[0]
    growth_err = d[d['parameter']=='max df std']['value'].values[0]
    dbl_time = d[d['parameter']=='inverse max df']['value'].values[0]
    dbl_err = d[d['parameter']=='inverse max df std']['value'].values[0]
   

    tidy_stats = tidy_stats.append({'date':g[0], 'carbon':g[1], 'temp_C':g[2], 'run_number':g[3],
                                    'growth_rate':growth_rate, 
                                    'dbl_time':dbl_time,
                                    'growth_err':growth_err,
                                    'dbl_err':dbl_err}, 
                                    ignore_index=True)
tidy_stats['growth_rate'] *= 60
tidy_stats['growth_err'] *= 60

#%%
# Compute the mean dbl time and growth rate over all replicates
growth = tidy_stats.groupby(['carbon', 'temp_C'])[
                            ['dbl_time', 'growth_rate']
                            ].agg(('mean', 'sem')).reset_index()



# %% Compute summary statistics of the fold-change data
data_grouped = data.groupby(['strain', 'carbon', 'temp', 'atc_ngml', 
                            'date', 'run_number']).mean().reset_index()
summarized = data_grouped.groupby(['strain', 'carbon', 'temp']).agg(('mean', 'sem')).reset_index()

#%%
fig, ax = plt.subplots(1, 1, dpi=100)
# ax.set_yscale('log')
delta = summarized[summarized['strain']=='delta']
for g, d in delta.groupby(['carbon', 'temp']):
    gr = growth[(growth['carbon']==g[0]) & (growth['temp_C']==g[1])]['growth_rate']['mean']
               
    ax.plot(gr, d['_sub']['mean'], '.', label=g)

ax.legend()
#%%
d

#%%
