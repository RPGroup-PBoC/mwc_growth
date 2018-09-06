import glob
import numpy as np
import pandas as pd
import os #maybe not needed. 
import sys
sys.path.insert(0, '../')
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
    if state['status'] == 'accepted':
        # Read data
        if expt[-1] != '/':
            expt += '/'
        _data = pd.read_csv('{}{}'.format(expt, 'output/*_flow_events.csv')

        # Record any extra or missing columns, and append data if no extra columns
        diff_extra = _data.columns.difference(columns) #any extra column names
        diff_miss = columns.difference(_data.columns) #any missing column names
        if (not diff_extra.empty) | (not diff_miss.empty):
            errors.append(expt.split('/')[-2] + 
                            ' -> extra: {}, missing: {} %\n'.format(
                            diff_extra.values.tolist(), diff_miss.values.tolist()))
        if diff_extra.empty:
            _dfs.append(_data)
            n_added += 1

df = pd.concat(_dfs, ignore_index=True)
df = df[columns.tolist()]

summary = """
Summary: ({} of {} experiments added)
=============================================================================
The following files are in error:\n
""".format(n_added, len(expt_paths))

if len(errors) == 0:
    failures = ['None\n']
for e in errors:
    summary += e

summary += """
-----------------------------------------------------------------------------
"""

# Save dataframe and error summary to output.
df.to_csv('../data/master_df.csv', index=False)
#what should the csv name be? include the date it was compiled or anything? 
#prob leave as is so new ones write over it? or maybe that's not ideal.

if os.path.isdir('../data/validation') == False: #where to put this summary? (we're sitting in code rn)
    os.mkdir('../data/validation')
with open('{}{}'.format('../data/validation/', 'flow_master_summary.txt'), 'w') as f:
    f.write(output)

print(summary) #what's proper here? nothing at all?