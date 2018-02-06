import numpy as np
import pandas as pd
import scipy.io
import tqdm


def clist_to_dataframe(clist_file, desired_props='default', added_props={},
                       excluded_props=None):
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
                         'Fluor2 mean death', 'Mother ID']
    defs = {key: value for value, key in enumerate(
        mat['def']) if key in desired_props}

    # Generate an empty DataFrame with the desired columns.
    for k, v in added_props.items():
        desired_props.append(k)
    df = pd.DataFrame([], columns=desired_props)

    # Iterate through the clist and extrac the properties.
    for i, cell in enumerate(mat['data']):
        # Extract the properties and add to DataFrame
        cell_dict = {key: cell[value] for key, value in defs.items()}

        # Add any additional properties.
        for k, v in added_props.items():
            cell_dict[k] = v
        df = df.append(cell_dict, ignore_index=True)

    # Rename the columns to accomodate pep8 style.
    if excluded_props is not None:
        for x in excluded_props:
            df.drop(x, axis=1, inplace=True)
    new_cols = {nom: '_'.join(nom.split(' ')).lower() for nom in df.keys()}
    df.rename(columns=new_cols, inplace=True)
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
        empty dict. If `parse_psotion` is True, the position will be passed as
        an added property.
    verbose: bool
        If True, a progressbar will be displayed for the clist iteration.

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
