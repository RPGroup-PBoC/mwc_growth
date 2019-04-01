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

PER_WELL = False

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
with open(f'output/{STRAIN}_{CARBON}_GP_summary.csv','w') as f:
    w = csv.writer(f)
    w.writerow(['parameter','value'])
    w.writerows(stats.items())
    
# Export full time series results of the fit
gp.export('output/gp_output.csv')

# Plot results
plt.figure(figsize=(10, 4))
plt.subplot(1,3,1)
plt.title('growth vs time', fontsize = 10)
gp.plotfit('f')
plt.subplot(1,3,2)
plt.title('time derivative vs time', fontsize=10)
gp.plotfit('df')
plt.subplot(1,3,3)
plt.title('doubling time vs time', fontsize=10)
gp.plotfit('dt')
plt.axhline(y=stats['inverse max df'],c='red')
locs, _ = plt.yticks()
plt.yticks(np.append(locs[1:-1],round(stats['inverse max df'],2)))
plt.tight_layout()
plt.savefig(f'output/{STRAIN}_{CARBON}_GP_fit_and_time_derivative.png')

# Perform by-well analysis
if PER_WELL:
    well_stats = []
    
    for well in restricted['well_id']:
        subset = restricted[restricted['well_id'] == well]
        
        # perform gaussian processing
        gp = mwc.fitderiv.fitderiv(subset['time_min'].values, subset['od_sub'].values)
        
        # export full results
        gp.export('output/well_by_well/gp_output_{well}.csv')

        # export statistics
        stats = gp.printstats(performprint=False)
        with open(f'output/per_well/{STRAIN}_{CARBON}_GP_summary_{well}.csv','w') as f:
            w = csv.writer(f)
            w.writerow(['parameter','value'])
            w.writerows(stats.items())
        
        # make dataframe of statistics and add to list
        _df = pd.DataFrame(stats)
        _df['well_id'] = well
        well_stats.append(_df)
        
        # plot and export results
        plt.figure(figsize=(10, 4))
        plt.subplot(1,3,1)
        plt.title('growth vs time', fontsize = 10)
        gp.plotfit('f')
        plt.subplot(1,3,2)
        plt.title('time derivative vs time', fontsize=10)
        gp.plotfit('df')
        plt.subplot(1,3,3)
        plt.title('doubling time vs time', fontsize=10)
        gp.plotfit('dt')
        plt.axhline(y=stats['inverse max df'],c='red')
        locs, _ = plt.yticks()
        plt.yticks(np.append(locs[1:-1],round(stats['doubling time'],2)))
        plt.tight_layout()
        plt.savefig(f'output/{STRAIN}_{CARBON}_GP_fit_and_time_derivative_{well}.png')
    
    # create single dataframe of the statistics for each well and save
    well_data = pd.concat(well_stats)
    well_data['carbon'] = CARBON
    well_data['strain'] = STRAIN
    well_data.to_csv('per_well_stats.csv', index=False)
    
# Troubleshooting
plt.figure(figsize=(10, 4))
plt.subplot(1,3,1)
plt.title('growth vs time', fontsize = 10)
gp.plotfit('f')
plt.subplot(1,3,2)
plt.title('time derivative vs time', fontsize=10)
gp.plotfit('df')
plt.subplot(1,3,3)
plt.title('doubling time vs time', fontsize=10)
gp.plotdoubtime()
plt.axhline(y=stats['doubling time'],c='red')
locs, _ = plt.yticks()
plt.yticks(np.append(locs[1:-1],round(stats['doubling time'],2)))
plt.tight_layout()
plt.savefig(f'output/{STRAIN}_{CARBON}_GP_fit_and_time_derivative_doubtime.png')