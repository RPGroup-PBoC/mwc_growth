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
expts = glob.glob('../../data/preprocessed/*')

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
    samps = glob.glob(f'{run}/snaps/*')
    for s in samps:
        strain, conc = s.split('/')[-1].split('_')
        conc = float(conc.split('ngml')[0])
        if strain is 'autoM':
            strain = 'auto' 

        clists = glob.glob(f'{s}/*/*.mat') 
        parsed = mwc.process.parse_clists(clists)

        # Assign the identifiers.
        parsed['strain'] = strain 
        parsed['carbon'] = carbon 
        parsed['atc_ngml'] = conc
        parsed['run_number'] = run_number
        parsed['temp'] = temp
        parsed['operator'] = operator 
        parsed['date'] = date
        intensities.append(parsed) 


# Concatenate all of the dataframes and save to disk. 
lineages = pd.concat(lineages)
intensities = pd.concat(intensities)
lineages.to_csv('../../data/raw_compiled_lineages.csv')
intensities.to_csv('../../data/raw_compiled_snaps.csv')
#%%
