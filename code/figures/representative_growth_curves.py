#%%
import numpy as np 
import pandas as pd 
import mwc.viz
import phd.viz  
import mwc.stats
import matplotlib.pyplot as plt
colors, color_list = mwc.viz.bokeh_theme()
_ = phd.viz.phd_style()
np.random.seed(666)

#%%
# Load the plates and statistics
plates = pd.read_csv('../../data/compiled_growth_plates.csv')
stats = pd.read_csv('../../data/compiled_growth_statistics.csv')
plates['temp_C'] = np.round(plates['temp_C'])
plates['time_min'] = np.round(plates['time_min'])
plates = plates[plates['time_min'] <= 600]

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

fig, ax = plt.subplots(1, 1, figsize=(3, 1.5), dpi=300)
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
ax.set_xticks([0, 100, 300, 500])

labels = ['acetate, 37°C', 'glycerol, 37°C','glucose, 37°C', 'glucose, 32°C', 'glucose, 42°C']
carbs = ['acetate', 'glycerol', 'glucose', 'glucose', 'glucose']
vals = ['acetate', 'glycerol', 'glucose', 32,  42]
temps = [37, 37, 37, 32, 42]
for c, t, l, v in zip(carbs, temps, labels, vals):
    for g, d in sub_plates[(sub_plates['carbon']==c) & (sub_plates['temp_C']==t)].groupby('carbon'):
        grp = d.groupby('time_min')['rel_od'].agg(('mean', 'std')).reset_index()
        ax.errorbar(grp['time_min'], grp['mean'], grp['std'], capsize=1.5, fmt='.', ms=3,
                lw=.75, label=l, linestyle='-', markerfacecolor=fill_colors[v],
                color=edge_colors[v], markeredgewidth=0.5)

ax.legend(handlelength=1, fontsize=8)
ax.set_xlabel('time [min]', style='italic', fontsize=10)
ax.set_ylabel('relative OD$_{600nm}$', style='italic', fontsize=10)
plt.savefig('../../figs/Fig1_growth_curve_comparison.pdf', bbox_inches='tight', facecolor='white')



#%%
# Instantiate the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(3, 1.5))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=8)

# Separate carbon and temperature
carbs = tidy_stats[(tidy_stats['temp_C']==37)].copy()
temps = tidy_stats[(tidy_stats['temp_C'] != 37)].copy()
carbs.loc[carbs['carbon']=='acetate', 'pos'] = 0
carbs.loc[carbs['carbon']=='glycerol', 'pos'] = 1
temps.loc[temps['temp_C']==42, 'pos'] = 3
temps.loc[temps['temp_C']==32, 'pos'] = 2

# Plot the scatters
for g, d in carbs.groupby(['carbon', 'pos']):
    if g[0] != 'glucose':
        ax.errorbar(np.random.normal(g[1], 0.1, len(d)), d['dbl_time'], d['dbl_err'], fmt='.', 
                alpha=0.75, markerfacecolor=fill_colors[g[0]], markeredgecolor=edge_colors[g[0]],
                markeredgewidth=0.75, ms=8, linestyle='none', capsize=1.5, color=edge_colors[g[0]])
for g, d in temps.groupby(['temp_C', 'pos']):
    ax.errorbar(np.random.normal(g[1], 0.1, len(d)), d['dbl_time'], d['dbl_err'], fmt='.', 
                alpha=0.75, markerfacecolor=fill_colors[g[0]], markeredgecolor=edge_colors[g[0]],
                markeredgewidth=0.75, ms=8, linestyle='none', capsize=1.5, color=edge_colors[g[0]])

glu = carbs[carbs['carbon']=='glucose']
ax.errorbar(np.random.normal(4, 0.1, len(glu)), glu['dbl_time'], glu['dbl_err'], fmt='.',
           alpha=0.75, markerfacecolor=fill_colors['glucose'], markeredgecolor=edge_colors['glucose'],
           color=edge_colors['glucose'], markeredgewidth=0.75, ms=8, linestyle='none', capsize=1.5)

ax.set_ylabel('doubling time [min]', style='italic', fontsize=10)
ax.set_xticks([0, 1, 2, 3, 4])
ax.set_xticklabels(['acetate\n37° C', 'glycerol\n37 °C', 'glucose\n37 °C',
                    'glucose\n42° C', 'glucose\n 32°C'], style='italic', rotation=45)
# ax.set_xlabel('growth condition', style='italic', fontsize=10)
ax.xaxis.grid(False)
ax.set_xlim([-0.5, 4.5])
ax.set_ylim([50, 225])
plt.savefig('../../figs/Fig1_growth_rate_jitter.pdf', facecolor='white', bbox_inches='tight')
#%%


#%%
