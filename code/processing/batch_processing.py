# ##############################################################################
# batchprocess.py
# ------------------------------------------------------------------------------
# Author: Griffin Chure
# License: MIT
#
# Description 
# ------------------------------------------------------------------------------
# This file does complete batch processing of all *approved* data sets for the
# project from pre-processed `clist.mat` files on local disk. It applies no area
# filters, no inference, no flattening -- just combines the data with relevant
# identifiers for further use.  
# ##############################################################################
#%%
import numpy as np
import pandas as pd
import glob
import mwc.process
import tqdm

# Instantiate storage lists for dataframes
lineages, intensities = [], []

# Define the list of experiments
expts = glob.glob('../../data/preprocessed/201*')

# Iterate through each experiment.
for run in tqdm.tqdm(expts, desc="Parsing experiments"):
    # Parse the info from the foldername
    date, run_number, temp, carbon, operator, _ = run.split('/')[-1].split('_')

    # Format the quantities as appropriate. 
    date = int(date)
    run_number = int(run_number.split('r')[-1])
    temp = int(temp.split('C')[0])

    # Parse all of the growth clists
    growth = glob.glob(f'{run}/growth/*/*.mat')

    # Reunite the families. 
    cells = mwc.process.parse_clists(growth)
    families = mwc.process.family_reunion(cells, fluo_channel=2)

    # Assign identifiers to the families
    families['carbon'] = carbon 
    families['run_number'] = run_number 
    families['temp'] = temp 
    families['date'] = date
    lineages.append(families)

    # Parse all of the snapshots
    samps = glob.glob(f'{run}/snaps/*ng*')
    _intensities = []
    for s in samps:
        strain, conc = s.split('/')[-1].split('_')
        conc = float(conc.split('ngml')[0])
        clists = glob.glob(f'{s}/*/*.mat') 
        parsed = mwc.process.parse_clists(clists)
        parsed = mwc.process.morphological_filter(parsed, ip_dist=0.065,
                                                area_bounds=[0.1, 4],
                                                ar_bounds=[0.1, 0.8])
        parsed['valid'] = parsed['error_frame'].isnull()
        parsed = parsed[parsed['valid']==True]
        if (strain == 'auto'):
            mean_auto_yfp = parsed['fluor1_mean_death'].mean() 
        if (strain == 'delta'):
            mean_delta_yfp = parsed['fluor1_mean_death'].values.mean()
        
        # Assign the identifiers.
        parsed['strain'] = strain 
        parsed['carbon'] = carbon 
        parsed['atc_ngml'] = conc
        parsed['run_number'] = run_number
        parsed['temp'] = temp
        parsed['operator'] = operator 
        parsed['date'] = date
        _intensities.append(parsed) 
    _intensities = pd.concat(_intensities)
    _intensities['fold_change'] = (_intensities['fluor1_mean_death'] - mean_auto_yfp) / \
                               (mean_delta_yfp - mean_auto_yfp) 
    intensities.append(_intensities)

# Concatenate all of the dataframes and save to disk. 
lineages = pd.concat(lineages)
intensities = pd.concat(intensities)

# Rename the columns to useful quantities
lineages = lineages[['I_1', 'I_2', 'area_1', 'area_2', 'parent_ID', 'sibling_ID_1', 'sibling_ID_2',
                        'carbon', 'run_number', 'temp', 'date', 'volume_1_birth', 
                        'volume_2_birth', 'volume_1_death', 'volume_2_death', 
                        'length_1_birth', 'length_2_birth',
                        'length_1_death', 'length_2_death', 'position']]
print(lineages.keys())
intensities = intensities[['area_death', 'fluor1_mean_death', 'fluor2_mean_death', 
                           'strain', 'date', 'run_number', 'temp',
                           'atc_ngml', 'carbon', 'long_axis_death',
                           'short_axis_death', 'volume_birth', 'volume_death', 
                           'area_birth',]]
intensities.rename(columns={'fluor1_mean_death': 'mean_yfp', 
                        'fluor2_mean_death':'mean_mCherry',
                        'area_death':'area_pix', 'long_axis_death':'length_um',
                        'short_axis_death':'width_um'}, inplace=True)
lineages.to_csv('../../data/raw_compiled_lineages.csv')
intensities.to_csv('../../data/raw_compiled_snaps.csv')
#%%
