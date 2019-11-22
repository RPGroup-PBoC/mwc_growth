import numpy as np
import pandas as pd
import scipy.io
import tqdm


def clist_to_dataframe(clist_file, desired_props='default', added_props={},
                       excluded_props=None, ip_dist=0.065):
    """
    Reads in a SuperSegger `clist` file and extracts the desired properties.

    Parameters
    ----------
    clist_file : str
        Path to clist file of interest.
    desired_props: list of str
        A list of the desired properties. Default selection is ['Area birth',
        'Area death', 'Cell ID', 'Cell birth time', 'Cell death time',
        'Daughter1 ID', 'Daughter2 ID', 'Fluor1 mean death',
        'Fluor2 mean death', 'Mother ID']
    added_props : dict
        A dict of additional props to include in the DataFrame.
    excluded_props : list or tuple of str
        Properties from the default included properties to ignore. These should
        match the case exactly as defined in the SuperSegger documentation.

    Returns
    -------
    df : pandas DataFrame
        A tidy pandas DataFrame with extracted properties for all cells in the
        clist file.
    """
    # Ensure that the clist file is a string.
    if type(clist_file) is not str:
        raise TypeError('clist_file must be a string')

    # Convert the excluded props to a list if is given as a string.
    if type(excluded_props) == str:
        excluded_props = list(excluded_props)

    # Load the clist file using scipy.
    mat = scipy.io.loadmat(clist_file, squeeze_me=True)

    # Assemble a dictionary of the indices and key values.
    if desired_props == 'default':
        desired_props = ['Area birth', 'Area death', 'Cell ID',
                         'Cell birth time', 'Cell death time', 'Daughter1 ID',
                         'Daughter2 ID', 'Fluor1 mean death',
                         'Fluor2 mean death', 'Mother ID', 
                         'Long axis (L) death', 'Long axis (L) birth',
                         'Short axis death', 'Short axis birth', 
                         'Fluor1 bg death', 'Fluor2 bg death',
                         'Error frame', 'Cell Dist to Edge']

    defs = {key: value for value, key in enumerate(
        mat['def']) if key in desired_props}
    # Generate an empty DataFrame with the desired columns.
    for k, v in added_props.items():
        desired_props.append(k)
    df = pd.DataFrame([], columns=desired_props)

    # Iterate through the clist and extrac the properties.
    for i, cell in enumerate(mat['data']):
        if type(cell) != np.float64:

            # Extract the properties and add to DataFrame
            cell_dict = {key: cell[value] for key, value in defs.items()}

            # Add any additional properties.
            for k, v in added_props.items():
                cell_dict[k] = v
            df = df.append(cell_dict, ignore_index=True)

    # Rename the columns to accommodate pep8 style.
    if excluded_props is not None:
        for x in excluded_props:
            df.drop(x, axis=1, inplace=True)
    new_cols = {nom: '_'.join(nom.split(' ')).lower() for nom in df.keys()}
    new_cols['Long axis (L) death'] = 'long_axis_death'
    new_cols['Long axis (L) birth'] = 'long_axis_birth'
    df.rename(columns=new_cols, inplace=True)

    # Compute the aspect ratio.
    df['short_axis_death'] = df['short_axis_death'] * ip_dist
    df['long_axis_death'] = df['long_axis_death'] * ip_dist
    df['short_axis_birth'] = df['short_axis_birth'] * ip_dist
    df['long_axis_birth'] = df['long_axis_birth'] * ip_dist

    df.loc[:, 'aspect_ratio'] = df['short_axis_death'] / df['long_axis_death']
    df.loc[:, 'volume_death'] = 0.5 * np.pi * df['short_axis_death'].values**2 *\
           ((2 * df['short_axis_death'].values / 3) + df['long_axis_death'].values
            - df['short_axis_death'].values)
    df.loc[:, 'volume_birth'] = 0.5 * np.pi * df['short_axis_birth'].values**2 *\
           ((2 * df['short_axis_birth'].values / 3) + df['long_axis_birth'].values
            - df['short_axis_birth'].values)

    return df

def parse_clists(clists, parse_position=True, added_props={},
                 verbose=False, **kwargs):
    """
    A helper function to iterate over a list of clist files. See
    `clist_to_dataframe` for function documentation.

    Parameters
    ----------
    clists: list of str
        A list of pathnames for clist files.
    parse_position: bool
        If True, the position of the item will be parsed from the file name by
        splitting at 'xy'.
    added_props: dict
        A dictionary of additional props to add to the DataFrame. Default is an
        empty dict. If `parse_position` is True, the position will be passed as
        an added property.
    verbose: bool
        If True, a progress bar will be displayed for the clist iteration.

    Returns
    -------
    df : pandas DataFrame
        A pandas DataFrame containing all cell properties for each item in the
        provided clist file.
    """

    # Iterate through each item in the clists.
    dfs = []
    if verbose:
        iterator = tqdm.tqdm(clists)
    else:
        iterator = clists
    for i, c in enumerate(iterator):
        # Parse the position.
        pos = int(c.split('xy')[-1].split('/')[0])
        if parse_position:
            added_props['position'] = pos

        # Pass the file to the parser.
        df = clist_to_dataframe(c, added_props=added_props, **kwargs)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def morphological_filter(df, ip_dist, area_bounds=[0, 50], ar_bounds=[0, 1]):
    """
    """
    # Calculate areas.
    area_bounds = np.array(area_bounds) / ip_dist**2

    _df = df.copy()
    _df = _df[(_df['area_death'] >= area_bounds[0]) &
              (_df['area_death'] <= area_bounds[1]) &
              (_df['aspect_ratio'] >= ar_bounds[0]) &
              (_df['aspect_ratio'] <= ar_bounds[1])]
    return _df


def family_reunion(dilution_df, multi_xy=True, fluo_channel=2):
    """"
    Generates a new DataFrame containing the individual sister intensities,
    summed fluorescence, and square fluctuations.

    Parameters
    ----------
    dilution_df : Pandas DataFrame
        The dataframe containing all lineage information. At a minimum,
        this must have the columns `mother_id`, `fluor1_mean_death`,
        `area_death`. If multiple positions are to be processed, `position`
        must also be provided.

    multi_xy : Bool
        If True, multiple positions in the provided DataFrame will be
        processed.

    Returns
    -------
    fluct_df : Pandas DataFrame
        A DataFrame with the two intensity measurements, sum total, and square fluctuations.
    """
    # Set up the DataFrame.
    family_df = pd.DataFrame([], columns=['I_1', 'I_2', 'area_1', 'area_2'
                                          'volume_1', 'volume_2'])

    # Determine what to groupby.
    if multi_xy == True:
        groupby = ['position', 'mother_id']
    else:
        groupby = ['mother_id']

    # Group the growth dataframe and iterate.
    grouped = dilution_df.groupby(groupby)
    for g, d in grouped:
        if multi_xy == True:
            mom = g[-1]
            pos = g[0]
        else:
            mom = g
            pos = 0
        if len(d) == 2:  # Ensure only single successful divisions.
            ints = d[f'fluor{fluo_channel}_mean_death'].values
            if sum(ints) > 0:
                I_1, I_2 = ints
                family_dict = {'I_1': I_1, 'I_2': I_2, 
                               'parent_ID': mom,
                               'sibling_ID_1': d['cell_id'].values[0],
                               'sibling_ID_2': d['cell_id'].values[1],
                               'error_frame': d['error_frame'].values[0],
                               'area_1':d['area_death'].values[0], 
                               'area_2':d['area_death'].values[1],
                               'position':pos,
                               'length_1_birth': d['long_axis_birth'].values[0],
                               'length_2_birth': d['long_axis_birth'].values[1],
                               'length_1_death': d['long_axis_death'].values[0],
                               'length_2_death': d['long_axis_death'].values[1],
                               'volume_1_birth': d['volume_birth'].values[0],
                               'volume_2_birth': d['volume_birth'].values[1],
                               'volume_1_death': d['volume_death'].values[0],
                               'volume_2_death': d['volume_death'].values[1],
                               'fractional_birth_area':d['area_birth'].values[0] / np.sum(d['area_birth'].values)}

                family_df = family_df.append(family_dict, ignore_index=True)

    return family_df

def condition_filter(df, strain='dilution', carbon=None, temp=None, fc_min=0, 
                  rep_min=0, size='large'):
    """
    Filters the fold-change data frames based on passed parameters. This
    function is used at the beginning of every script to ensure all data is
    being processed the same way. 

    Parameters
    ----------
    df : pandas DataFrame
        The dataframe to be filtered.
    strain : str, NoneType
        The desired strain to be isolated. Default is 'dilution'. If None is
        passed all strains will be returned. 
    carbon: str, NoneType
        The desired carbon source to be isolated. Default is None
        and all carbon sources will be returned. 
    temp: int, NoneType
        The desired temperature to be isolated. Default is None
        and all temperatures will be returned. 
    fc_min: float   
        The minimum (noise floor, inclusive) fold-change to be considered. Default is 0 
    rep_min: float
        The minimum (noise floor, exclusive) number of repressors to be considered. 
        Default is 0.
    size: str, ['large', 'medium', 'small']
        Desired size of cells to be returned. Default is "large"

    Returns
    -------
    filt_df : pandas DataFrame
        A copy of the dataframe with filters applied. 
    """

    if strain not in [None, 'dilution', 'delta', 'auto']:
       raise ValueError(f"Provided strain is {strain}. Must be either 'dilution', 'auto', 'delta', or None")
    if carbon not in [None, 'glucose', 'glycerol', 'acetate']:
       raise ValueError(f"Provided carbon is {carbon}. Must be either 'glucose', 'glycerol', 'acetate', or None")
    if temp not in [None, 37, 32, 42]:
       raise ValueError(f"Provided carbon is {temp}. Must be either 32, 37, 42, or None")
    
    # Make a copy of the data frame and apply fold-chagne and repressor filters.
    filt_df = df.copy()
    filt_df = filt_df[(filt_df['fold_change'] >= fc_min) & 
              (filt_df['repressors'] > rep_min)]
    
    # Determine if the other filters should be applied. 
    if carbon != None:
        filt_df = filt_df[filt_df['carbon']==carbon]
    if temp != None:
        filt_df = filt_df[filt_df['temp']==temp]
    if strain != None:
        filt_df = filt_df[filt_df['strain']==strain]
    if size != None:
        filt_df = filt_df[filt_df['size']==size]
    
    return filt_df