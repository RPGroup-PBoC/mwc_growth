import glob
import numpy as np
import pandas as pd
import os #maybe not needed. 
import sys
sys.path.insert(0, '../')
sys.path.insert(0, '../../work/PhillipsLab/git/mwc_growth')
import mwc.io

# Get the experiment paths.
expt_paths = glob.glob('./processing/*_dilution/')
expt_paths = np.sort(expt_paths)


columns = pd.Index(['date', 'run', 'operator', 'carbon', 'temp', 'strain',
                    'atc_ngml', 'iptg_um', 'mean_fitc','fold_change','reps'])
_dfs = []
errors = []
n_added = 0

# Append data to list, then concat into dataframe.
for i, expt in enumerate(expt_paths):
    # Check status 
    state = mwc.io.scrape_frontmatter(expt)
    if len(state) == 0:
        print('skipped')
    elif state['status'].lower() == 'accepted':
        # Read data
        if expt[-1] != '/':
            expt += '/'
        for f in glob.glob('{}output/*_flow_events.csv'.format(expt)): #expect only 1 file
            _data = pd.read_csv(f)

        # Record any extra or missing columns, and append data if no extra columns
        diff_extra = _data.columns.difference(columns) #any extra column names
        diff_miss = columns.difference(_data.columns) #any missing column names
        if (not diff_extra.empty) | (not diff_miss.empty):
            errors.append('    ' + expt.split('/')[-2] + 
                            ' -> extra: {}, missing: {} \n'.format(
                            diff_extra.values.tolist(), diff_miss.values.tolist()))
        if diff_extra.empty:
            _dfs.append(_data)
            n_added += 1

if (len(_dfs) != 0):
    df = pd.concat(_dfs, ignore_index=True)
    if (set(df.columns) == set(columns)):
        df = df[columns.tolist()]
    else:
        print('Column(s) are missing, so column ordering was skipped.')

    # Create error summary
    summary = """
    Summary: ({} of {} experiments added)
    =============================================================================
    The following files are in error:\n \n""".format(n_added, len(expt_paths))

    if len(errors) == 0:
        failures = ['None\n']
    for e in errors:
        summary += e

    summary += \
    """    -----------------------------------------------------------------------------
    """

    # Save dataframe and error summary to output.
    df.to_csv('../data/master_df.csv', index=False)

    if os.path.isdir('../data/validation') == False:
        os.mkdir('../data/validation')
    with open('{}{}'.format('../data/validation/', 'flow_master_summary.txt'), 'w') as f:
        f.write(summary)

    print(summary)

else:
    print('No dataframes were compiled.')