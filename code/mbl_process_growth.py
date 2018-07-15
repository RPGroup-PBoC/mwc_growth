# %%

import numpy as np
import pandas as pd
import sys
sys.path.insert(0, '../')
import glob
import mwc.process
DATE = 20180621
RUN = 1
TEMP = 37
CARBON = 'glucose'
OPERATOR = 'O2'

BASE_NAME = '{}_{}_{}C_{}_{}'.format(DATE, RUN, TEMP, CARBON, OPERATOR)

IP_DIST = 0.065 # interpixel distance
# %% Processing of data
# data_dir = '/Volumes/{}_dilution/'.format(BASE_NAME)

# Extract growth (timelapse) file names and parse.
growth_files = glob.glob('/Volumes/GDC_DATA_2/20180621_growth/xy*/clist.mat')
print(growth_files)
growth_df = mwc.process.parse_clists(growth_files)


# Apply area bounds filter.
growth_df['no_err'] = growth_df['error_frame'].isnull() # returns a bool (1 or 0)
growth_df = growth_df[growth_df['no_err'] == 1]
growth_df = growth_df[growth_df['fluor1_mean_death'] > 0]


def compute_fluctuations(dilution_df, multi_xy=True,
                         fluo_key='fluor1_mean_death'):
    """"
    Generates a new DataFrame containing the indivudual sister intensities,
    summed fluorescence, and square fluctuations.

    Parameters
    ----------
    dilution_df : Pandas DataFrame
        The dataframe containing all lineage information. At a minimum,
        this must have the columns `mother_id`, `fluor1_mean_death`,
        `area_death`. If multiple positions are to be processed, `position`
        must also be provided.
    auto_val : float
        The mean autofluorescence pixel intensity.
    multi_xy : Bool
        If True, multiple positions in the provided DataFrame will be
        processed.

    Returns
    -------
    fluct_df : Pandas DataFrame
        A DataFrame with the two intensity measurements, sum total, and square fluctuations.
    """
    # Set up the DataFrame.
    fluct_df = pd.DataFrame([], columns=['I_1', 'I_2'])

    # Determine what to groupby.
    if multi_xy == True:
        print('Grouping by position and mother ID....')
        groupby = ['position', 'mother_id']
    else:
        print('Grouping only by mother ID...')
        groupby = ['mother_id']

    # Group the growth dataframe and iterate.
    grouped = dilution_df.groupby(groupby)
    for g, d in grouped:
        if len(d) == 2:  # Ensure only single successful divisions.
            ints = d[fluo_key].values * d['area_death'].values
            if (ints >= 0).all() == True:
                I_1, I_2 = ints
                summed = np.sum(ints)
                fluct = (ints[0] - ints[1])**2
                family_dict = {'I_1': I_1, 'I_2': I_2}

                fluct_df = fluct_df.append(family_dict, ignore_index=True)

    return fluct_df

fluct_df = compute_fluctuations(growth_df)
fluct_df.to_csv('{}.csv'.format(BASE_NAME), index=None)