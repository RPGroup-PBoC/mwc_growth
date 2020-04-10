#%%
import numpy as np 
import pandas as pd 
import mwc.viz
import mwc.stats
import matplotlib.pyplot as plt
colors, color_list = mwc.viz.plotting_style()

#%%
# Load the plates and statistics
plates = pd.read_csv('../../data/compiled_growth_plates.csv', comment='#')
stats = pd.read_csv('../../data/compiled_growth_statistics.csv', comment="#")
plates['temp_C'] = np.round(plates['temp_C'])
plates['time_min'] = np.round(plates['time_min'])
plates = plates[(plates['time_min'] <= 600)]

#%%
# reform the stats 
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
# Restrict only to conditions studied here.
plates = plates[((plates['carbon']=='LB') | (plates['carbon']=='glucose') | (plates['carbon']=='acetate') |
                (plates['carbon']=='glycerol') | (plates['carbon']=='blank')) & 
                ((plates['temp_C']==37) |  (plates['temp_C']==32) | 
                (plates['temp_C']==42))] 

#%% Subtract the blank from each date, run, and time point.
sub_plates = []
for g, d in plates.groupby(['date', 'run_number', 'time_min']):
    d = d.copy()
    blank = d[d['carbon']=='blank']['od_600nm'].mean()
    d['od_sub'] = d['od_600nm'] - blank
    sub_plates.append(d)
sub_plates = pd.concat(sub_plates)
_sub_plates = []
for g, d in sub_plates.groupby(['date', 'carbon', 'temp_C']):
    d = d.copy()
    d = d[d['od_sub']>0]
    a0 = d[d['time_min']==6]['od_sub'].mean()
    d['rel_od'] = d['od_sub'] / a0
    _sub_plates.append(d)
sub_plates = pd.concat(_sub_plates)
#%%
# Choose representable samples. 
#Define the colors
fill_colors = {'acetate': '#e1bb96', 'glycerol': colors['light_green'],
               'glucose':colors['light_purple'], 37: colors['light_purple'],
               32:colors['light_blue'], 42:colors['light_red']}
edge_colors = {'acetate': '#764f2a', 'glycerol': colors['green'],
               'glucose':colors['purple'], 37: colors['purple'],
               32:colors['blue'], 42:colors['red']}

fig, ax = plt.subplots(1, 1, figsize=(3.25, 3), dpi=300)
mwc.viz.despine(ax)
ax.set_yscale('log')
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
ax.set_xticks([0, 100, 300, 500])
ax.set_yticks([1, 10, 100])
ax.set_xlim([100, 500])
ax.set_ylim([1, 100])

labels = ['acetate, 37°C', 'glycerol, 37°C','glucose, 37°C', 'glucose, 32°C', 'glucose, 42°C']
carbs = ['acetate', 'glycerol', 'glucose', 'glucose', 'glucose']
vals = ['acetate', 'glycerol', 'glucose', 32,  42]
temps = [37, 37, 37, 32, 42]
for c, t, l, v in zip(carbs, temps, labels, vals):
    for g, d in sub_plates[(sub_plates['carbon']==c) & (sub_plates['temp_C']==t)].groupby('carbon'):
        grp = d.groupby('time_min')['rel_od'].agg(('mean', 'std')).reset_index()
        ax.errorbar(grp['time_min'], grp['mean'], grp['std'], capsize=1, fmt='o', ms=3,
               lw=0.5, label='__nolegend__', linestyle='', markeredgecolor='white',
                color=edge_colors[v], markeredgewidth=0.25)
for c, t, l, v in zip(carbs, temps, labels, vals):
    _stats = stats[(stats['carbon']==c) & (stats['temp']==t) & 
                   (stats['parameter']=='max df')]['value'].agg(('mean', 'sem'))
    dbl_mean = np.log(2) / _stats['mean']
    dbl_max = (np.log(2) / (_stats['mean'] - _stats['sem'])) - dbl_mean
    t_max = sub_plates[(sub_plates['carbon']==c) & 
                            (sub_plates['temp_C']==t) & 
                            (sub_plates['time_min']==503)]['rel_od'].mean()
    if c == 'glycerol':
        fudge = - 0.1
    else:
        fudge = 0 
    ax.text(505, t_max - (0.15 - fudge) * t_max, 
        r"$t_\mathrm{double} = $" + 
        f"{int(np.round(dbl_mean, decimals=0))} $\pm$ {int(np.round(dbl_max, decimals=0))} min",
        fontsize=6, color=edge_colors[v])
    ax.text(505, t_max + (0.05 + fudge) * t_max, l, fontsize=6, color=edge_colors[v])

# ax.legend(handlelength=1, fontsize=6)
ax.set_xlabel('time [min]', fontsize=8)
ax.set_ylabel('relative OD$_{600nm}$', fontsize=8)
plt.savefig('../../figs/Fig1_growth_curves.pdf', facecolor='none')



# %%


# %%
