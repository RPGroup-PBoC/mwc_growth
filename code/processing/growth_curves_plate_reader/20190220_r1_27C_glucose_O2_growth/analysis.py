# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../../../')
import numpy as np
import pandas as pd
import string
import matplotlib.pyplot as plt
import mwc.viz
colors = mwc.viz.personal_style()
import mwc.fitderiv
import csv

# Define the experimental constants
STRAIN = 'delta'
CARBON = 'glucose'

PER_WELL = True

# ----------------------------------

# Load the data. 
data = pd.read_csv(f'output/growth_plate.csv')

# Generate a dictionary of the mean blank at each time point. 
blank_vals = {t:val['od_600nm'].mean() for t, val in data[data['strain']=='blank'].groupby(['time_min'])}

# Add mean blank values for each time point to the dataframe, as well as background subtracted OD values.
for k, v in blank_vals.items():
    data.loc[data['time_min']==k, 'blank_val'] = v
data['od_sub'] = data['od_600nm'] - data['blank_val']

# Load time range to be analyzed
ranges_df = pd.read_csv('output/ranges.csv')
XMIN = ranges_df[STRAIN+CARBON][0]
XMAX = ranges_df[STRAIN+CARBON][1]

# Restrict data to desired range
restricted = data[(data['strain']==STRAIN) & (data['carbon']==CARBON) 
                  & (data['time_min'] > XMIN) & (data['time_min'] < XMAX)].sort_values('time_min')

# Using the package [`fitderiv`](http://swainlab.bio.ed.ac.uk/software/fitderiv/) from Peter Swain's lab, 
# perform non-parametric inference of the time-dependent growth rates.
gp = mwc.fitderiv.fitderiv(restricted['time_min'].values, restricted['od_sub'].values)

# Export summary statistics
stats = gp.printstats(performprint=False)
with open(f'output/{STRAIN}_{CARBON}/gp_output_stats.csv','w') as f:
    w = csv.writer(f)
    w.writerow(['parameter','value'])
    w.writerows(stats.items())
    
# Export full time series results of the fit
gp.export(f'output/{STRAIN}_{CARBON}/gp_output.csv')

# Read in gp results as dataframe and calculate doubling time
gp_df = pd.read_csv(f'output/{STRAIN}_{CARBON}/gp_output.csv')
gp_df = gp_df[['t','od','log(OD)','log(OD) error','gr','gr error']]
gp_df.rename(columns = {'log(OD)':'log(OD)_fit', 'log(OD) error':'log(OD)_fit_std', 
                        'gr':'growth_rate', 'gr error':'growth_rate_std', 
                        'od':'OD_raw_data', 't':'time'}, inplace = True)
gp_df['doubling_time'] = np.log(2)/gp_df['growth_rate']
gp_df['doubling_time_std'] = np.log(2)*gp_df['growth_rate_std']/(gp_df['growth_rate']**2)
gp_df.to_csv(f'output/{STRAIN}_{CARBON}/gp_output.csv')

# Plot results
fig, ax = plt.subplots(ncols = 3, figsize=(10, 4))
ax[0].set_title('growth vs time', fontsize = 10)
ax[1].set_title('time derivative vs time', fontsize = 10)
ax[2].set_title('doubling time vs time', fontsize = 10)

ax[0].scatter(gp_df['time'],np.log(gp_df['OD_raw_data']),c='r', marker = '.')
ax[0].plot(gp_df['time'],gp_df['log(OD)_fit'],c='blue')
ax[0].fill_between(gp_df['time'], gp_df['log(OD)_fit']-gp_df['log(OD)_fit_std'], 
                      gp_df['log(OD)_fit']+gp_df['log(OD)_fit_std'],
                      facecolor= 'blue', alpha= 0.2)

ax[1].plot(gp_df['time'],gp_df['growth_rate'],c='b')
ax[1].fill_between(gp_df['time'], gp_df['growth_rate']-gp_df['growth_rate_std'], 
                      gp_df['growth_rate']+gp_df['growth_rate_std'],
                      facecolor= 'blue', alpha= 0.2)

ax[2].plot(gp_df['time'],gp_df['doubling_time'],c='b')
ax[2].fill_between(gp_df['time'], gp_df['doubling_time']-gp_df['doubling_time_std'], 
                      gp_df['doubling_time']+gp_df['doubling_time_std'],
                      facecolor= 'blue', alpha= 0.2)
ax[2].axhline(y=gp_df.min()['doubling_time'],c='red')
locs = ax[2].get_yticks()
plt.yticks(np.append(locs[2:-1],round(gp_df.min()['doubling_time'],2)))

plt.tight_layout()
plt.savefig(f'output/{STRAIN}_{CARBON}/gp_output_curves.png')

# Perform by-well analysis
if PER_WELL:
    well_stats = []
    
    for well in restricted['well_id'].unique():
        # select data for single well
        subset = restricted[restricted['well_id'] == well]
        
        # run gaussian processing and export fit and time derivative results
        gp = mwc.fitderiv.fitderiv(subset['time_min'].values, subset['od_sub'].values)
        gp.export(f'output/{STRAIN}_{CARBON}/per_well/gp_output_{well}.csv')
        
        # export summary statistics
        stats = gp.printstats(performprint=False)
        with open(f'output/{STRAIN}_{CARBON}/per_well/gp_output_stats_{well}.csv','w') as f:
            w = csv.writer(f)
            w.writerow(['parameter','value'])
            w.writerows(stats.items())
        
        # add summary statistics to list
        stats['well_id'] = well
        well_stats.append(pd.DataFrame(stats))
        
        # read in gp results exported above as dataframe and calculate doubling time
        gp_df = pd.read_csv(f'output/{STRAIN}_{CARBON}/per_well/gp_output_{well}.csv')
        gp_df = gp_df[['t','od','log(OD)','log(OD) error','gr','gr error']]
        gp_df.rename(columns = {'log(OD)':'log(OD)_fit', 'log(OD) error':'log(OD)_fit_std', 
                                'gr':'growth_rate', 'gr error':'growth_rate_std', 
                                'od':'OD_raw_data', 't':'time'}, inplace = True)
        gp_df['doubling_time'] = np.log(2)/gp_df['growth_rate']
        gp_df['doubling_time_std'] = np.log(2)*gp_df['growth_rate_std']/(gp_df['growth_rate']**2)
        gp_df.to_csv(f'output/{STRAIN}_{CARBON}/per_well/gp_output.csv')
        
        # plot results
        fig, ax = plt.subplots(ncols = 3, figsize=(10, 4))
        ax[0].set_title('growth vs time', fontsize = 10)
        ax[1].set_title('time derivative vs time', fontsize = 10)
        ax[2].set_title('doubling time vs time', fontsize = 10)

        ax[0].scatter(gp_df['time'],np.log(gp_df['OD_raw_data']),c='r', marker = '.')
        ax[0].plot(gp_df['time'],gp_df['log(OD)_fit'],c='blue')
        ax[0].fill_between(gp_df['time'], gp_df['log(OD)_fit']-gp_df['log(OD)_fit_std'], 
                              gp_df['log(OD)_fit']+gp_df['log(OD)_fit_std'],
                              facecolor= 'blue', alpha= 0.2)

        ax[1].plot(gp_df['time'],gp_df['growth_rate'],c='b')
        ax[1].fill_between(gp_df['time'], gp_df['growth_rate']-gp_df['growth_rate_std'], 
                              gp_df['growth_rate']+gp_df['growth_rate_std'],
                              facecolor= 'blue', alpha= 0.2)

        ax[2].plot(gp_df['time'],gp_df['doubling_time'],c='b')
        ax[2].fill_between(gp_df['time'], gp_df['doubling_time']-gp_df['doubling_time_std'], 
                              gp_df['doubling_time']+gp_df['doubling_time_std'],
                              facecolor= 'blue', alpha= 0.2)
        ax[2].axhline(y=gp_df.min()['doubling_time'],c='red')
        locs = ax[2].get_yticks()
        plt.yticks(np.append(locs[2:-1],round(gp_df.min()['doubling_time'],2)))
        plt.tight_layout()
        plt.savefig(f'output/{STRAIN}_{CARBON}/per_well/gp_output_curves_{well}.png')
    
    # compile dataframe of statistics for all wells
    well_data = pd.concat(well_stats)
    well_data['carbon'] = CARBON
    well_data['strain'] = STRAIN
    well_data.to_csv(f'output/{STRAIN}_{CARBON}/per_well/per_well_stats.csv', index=False)