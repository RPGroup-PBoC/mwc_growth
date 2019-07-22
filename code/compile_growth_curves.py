#%%
import numpy as np
import pandas as pd 
import os 
import glob
import mwc.io
import tqdm

# Get a list of readme files 
experiments = glob.glob('./processing/growth_curves/2019*') 

# Isolate the accepted experiments. 
accepted = []
for exp in experiments:
    info = mwc.io.scrape_frontmatter(exp)
    if info['status'].lower() == 'accepted':
        accepted.append(exp)

# Instantiate the storage vectors 
gp_stats = []
plates = []

# Iterate through each directory and go to town. 
for i, exp in enumerate(tqdm.tqdm(accepted)):
    # Parse date and run information
    date, run_number, temp, _, operator, _ = exp.split('/')[-1].split('_')
    run_number = int(run_number.split('r')[-1])
    temp = float(temp.split('C')[0])
    date = int(date)

    # Append the growth plate file 
    plates.append(pd.read_csv(f'{exp}/output/growth_plate.csv'))
    
    # Get a list of the 'delta_' files which are the experimental strains. 
    samps = glob.glob(f'{exp}/output/delta_*')
    for s in samps:
        # Parse the carbon  and file
        carbon = s.split('/')[-1].split('_')[-1]
        output = pd.read_csv(f'{s}/gp_output_stats.csv')

        # Add identifying information
        output['date'] = date
        output['run_number'] = run_number
        output['temp'] = temp 
        output['carbon'] = carbon

        # Append to the list
        gp_stats.append(output)

# Concatenate and drop unnecessary columns
plates = pd.concat(plates)
stats = pd.concat(gp_stats)
plates = plates[['carbon', 'date', 'run_number', 'strain', 
                        'temp_C', 'time_min', 'well_id',
                        'od_600nm']]
stats = stats[['parameter', 'value', 'date', 
               'run_number','temp', 'carbon']]

# Keep only the true experimental strains and the blank
plates = plates[(plates['strain']=='delta') | (plates['strain']=='blank')]
plates.to_csv('../data/compiled_growth_plates.csv', index=False)
stats.to_csv('../data/compiled_growth_statistics.csv', index=False)


#%%


#%%
