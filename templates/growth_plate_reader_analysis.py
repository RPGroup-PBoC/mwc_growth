# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../../../')
import numpy as np
import pandas as pd
import string
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mwc.viz
colors = mwc.viz.personal_style()
import mwc.fitderiv
import seaborn as sns
import statsmodels.api as sm

# Define the experimental constants
# For figure titles
DATE = 
RUN_NO = 
TEMP = 
# For analysis
STRAIN = ''
CARBON = ''

WHOLE_PLATE = True
WHOLE_REPLOT_ONLY = False #only replots if whole plate is also true
PER_WELL = True
PER_WELL_REPLOT_ONLY = False

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

# Make output directory
if os.path.exists(f'output/{STRAIN}_{CARBON}/') == False:
    os.mkdir(f'output/{STRAIN}_{CARBON}/')

# Perform full plate analysis
if WHOLE_PLATE:
    if not WHOLE_REPLOT_ONLY:
        # Using the package [`fitderiv`](http://swainlab.bio.ed.ac.uk/software/fitderiv/) from Peter Swain's lab, 
        # perform non-parametric inference of the time-dependent growth rates.
        gp = mwc.fitderiv.fitderiv(restricted['time_min'].values, restricted['od_sub'].values)

        # Export summary statistics
        stats = gp.printstats(performprint=False)
        _stats = pd.DataFrame(list(stats.items())).rename(columns = {0:'parameter',1:'value'})
        _stats.to_csv(f'output/{STRAIN}_{CARBON}/gp_output_stats.csv')

        # Create dataframe with full time series results of the fit
        gp_df = gp.export('NONE', savegp = False, savestats = False)
        # Calculate doubling time
        gp_df = gp_df[['t','od','log(OD)','log(OD) error','gr','gr error']]
        gp_df.rename(columns = {'log(OD)':'log(OD)_fit', 'log(OD) error':'log(OD)_fit_std', 
                                'gr':'growth_rate', 'gr error':'growth_rate_std', 
                                'od':'OD_raw_data', 't':'time'}, inplace = True)
        gp_df['doubling_time'] = np.log(2)/gp_df['growth_rate']
        gp_df['doubling_time_std'] = np.log(2)*gp_df['growth_rate_std']/(gp_df['growth_rate']**2)
        # Save results
        gp_df.to_csv(f'output/{STRAIN}_{CARBON}/gp_output.csv')
    else:
        gp_df = pd.read_csv(f'output/{STRAIN}_{CARBON}/gp_output.csv')
    
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
    if not PER_WELL_REPLOT_ONLY:
        well_stats = []
        well_gps = []

        for well in restricted['well_id'].unique():
            # Select data for single well
            subset = restricted[restricted['well_id'] == well]

            # Run gaussian processing and export fit and time derivative results
            gp = mwc.fitderiv.fitderiv(subset['time_min'].values, subset['od_sub'].values)
            gp_df = gp.export('NONE', savegp = False, savestats = False)

            # Append summary statistics
            stats = gp.printstats(performprint=False)
            well_stats.append(pd.DataFrame(list(stats.items())).assign(well_id = well))

            # Calculate doubling time and add to gp results, then append dataframe to list of gp results
            gp_df = gp_df.assign(well_id = well)
            gp_df = gp_df[['well_id','t','od','log(OD)','log(OD) error','gr','gr error']]
            gp_df.rename(columns = {'log(OD)':'log(OD)_fit', 'log(OD) error':'log(OD)_fit_std', 
                                    'gr':'growth_rate', 'gr error':'growth_rate_std', 
                                    'od':'OD_raw_data', 't':'time'}, inplace = True)
            gp_df['doubling_time'] = np.log(2)/gp_df['growth_rate']
            gp_df['doubling_time_std'] = np.log(2)*gp_df['growth_rate_std']/(gp_df['growth_rate']**2)
            well_gps.append(gp_df)
        
        # Compile dataframe of GP results for all wells
        well_data = pd.concat(well_gps)
        well_data['carbon'] = CARBON
        well_data['strain'] = STRAIN
        well_data = well_data.sort_values('well_id')
        well_data.to_csv(f'output/{STRAIN}_{CARBON}/per_well_output.csv', index=False)

        # Compile dataframe of statistics for all wells
        well_stats_df = pd.concat(well_stats)
        well_stats_df = well_stats_df.rename(columns = {0:'parameter',1:'value'})
        well_stats_df['carbon'] = CARBON
        well_stats_df['strain'] = STRAIN
        well_stats_df = well_stats_df.sort_values('well_id')
        well_stats_df.to_csv(f'output/{STRAIN}_{CARBON}/per_well_stats.csv', index=False)

    else:
        well_data = pd.read_csv(f'output/{STRAIN}_{CARBON}/per_well_output.csv')
        well_stats_df = pd.read_csv(f'output/{STRAIN}_{CARBON}/per_well_stats.csv')

    alpha_map = {alpha:no for alpha, no in zip(string.ascii_uppercase, np.arange(0, 27, 1) * 12)}
    alphanumeric_map = {f'{a}{n}':alpha_map[a] + n for n in np.arange(1, 13, 1) for a in string.ascii_uppercase}
    row_letters = {no:alpha for no, alpha in zip(np.arange(0, 27, 1),string.ascii_uppercase)}

    # Plot doubling time curves for all wells
    fig, ax = plt.subplots(8,12, figsize=(10, 4), sharex=True) #,sharey=True

    for r in np.arange(0,8,1):
        for c in np.arange(0,12,1):
            well_id = alphanumeric_map[f'{row_letters[r]}{c+1}']
            data = well_data[well_data['well_id']==well_id].sort_values('time')
            ax[r][c].get_xaxis().set_visible(False)
            ax[r][c].get_yaxis().set_visible(False)
            if not data.empty:
                ax[r][c].plot(data['time'],data['doubling_time'],c='b')
                ax[r][c].fill_between(data['time'], data['doubling_time']-data['doubling_time_std'], 
                                      data['doubling_time']+data['doubling_time_std'],
                                      facecolor= 'blue', alpha= 0.2)
                ax[r][c].axvline(x=int(data[data['doubling_time'] == data['doubling_time'].min()]['time']),c='red')
                ax[r][c].set_ylim([0, 300])
    ax[0][0].set_ylim([0, 300])
    ax[0][0].get_yaxis().set_visible(True)
    fig.suptitle(f'{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{STRAIN} doubling time vs time', fontsize = 10, y=0.93)
    plt.savefig(f'output/{STRAIN}_{CARBON}/per_well_doubling_time_curves.png')    
    plt.close()

    # Plot heatmap of minimum doubling times for all wells
    fig, ax = plt.subplots()
    plate = np.full((8, 12),np.nan)
    for r in np.arange(0,8,1):
        for c in np.arange(0,12,1):
            well_id = alphanumeric_map[f'{row_letters[r]}{c+1}']
            _temp = well_stats_df[well_stats_df['well_id']==well_id]
            if not _temp.empty:
                plate[r][c] = int(_temp[_temp['parameter'] == 'inverse max df']['value'])
    mask = np.isnan(plate)
    ax = sns.heatmap(plate, mask=mask, square = True)
    ax.set_title(f'{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{STRAIN} doubling times', fontsize = 10)
    plt.savefig(f'output/{STRAIN}_{CARBON}/per_well_doubling_times_heatmap.png')    
    plt.close()

    # Plot distribution of minimum doubling times for all wells
    fig, ax = plt.subplots()
    dat = well_stats_df[well_stats_df['parameter']=='inverse max df']['value']
    ecdf = sm.distributions.ECDF(dat)
    x = np.linspace(min(dat), max(dat))
    y = ecdf(x)
    ax.step(x, y)
    ax.set_xlabel('doubling time (minutes)', fontsize = 10)
    ax.set_ylabel('frequency', fontsize = 10)
    ax.set_title(f'{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{STRAIN} doubling time distribution across wells', fontsize = 10)
    ax.tick_params(labelsize=10)
    plt.savefig(f'output/{STRAIN}_{CARBON}/per_well_doubling_times_distribution.png')
    plt.close()