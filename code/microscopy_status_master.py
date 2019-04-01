import glob
import numpy as np
import pandas as pd
import os
import sys
sys.path.insert(0, '../')
import mwc.io

# Get the experiment paths.
expt_paths = glob.glob('./processing/microscopy/*_dilution/')
expt_paths = np.sort(expt_paths)

basename, status, reason, errors, accepted, rejected, no_status, carbon = \
    [ [] for i in range(8)]
n_added, n_errors = [0,0]

# Append names and states to respective lists, then concat into dataframe.
for i, expt in enumerate(expt_paths):
    name = expt.split('/')[-2]

    # Check status 
    state = mwc.io.scrape_frontmatter(expt)
    if len(state) == 0:
        print('skipped')
        no_status.append(name)
    elif state['status'].lower() == 'accepted':
        accepted.append(name)
    elif state['status'].lower() == 'rejected':
        rejected.append(name)
    else:
        n_errors += 1
        errors.append(name)
    n_added += 1
    basename.append(name)
    carbon.append(name.split('_')[3])
    status.append(state['status'])
    reason.append(state['reason'])

if (len(basename) != 0):
    df = pd.DataFrame({'experiment':basename, 'carbon':carbon, 'status':status, 'reason':reason})
    df = df[['experiment','carbon','status','reason']]
    df.sort_values('carbon').sort_values('reason')

    # Create summary
    summary = """
    Summary: ({} of {} experiments added)
    =============================================================================
    The following files {} are in error:\n \n""".format(n_added, len(expt_paths), n_errors)

    if len(errors) == 0:
        summary += '    None\n'
    for i in errors:
        summary += ('    '+i+'\n')

    summary += \
    """    -----------------------------------------------------------------------------
    The following files {} were accepted:\n \n""".format(len(accepted))

    if len(accepted) == 0:
        summary += '    None\n'
    summary += '    Glucose:\n'
    for i in accepted:
        carbon_val = i.split('_')[3]
        date = int(i.split('_')[0])
        if carbon_val == 'glucose':
            summary += ('    '+i)
            if date > 20181126: summary += (' (.5%)\n')
            else: summary += ('\n')
    summary += '    Glycerol:\n'
    for i in accepted:
        carbon_val = i.split('_')[3]
        date = int(i.split('_')[0])
        if carbon_val == 'glycerol':
            summary += ('    '+i)
            if date > 20181126: summary += (' (.5%)\n')
            else: summary += ('\n')
    summary += '    Acetate:\n'
    for i in accepted:
        carbon_val = i.split('_')[3]
        date = int(i.split('_')[0])
        if carbon_val == 'acetate':
            summary += ('    '+i)
            if date > 20181126: summary += (' (.5%)\n')
            else: summary += ('\n')

    summary += \
    """    -----------------------------------------------------------------------------
    The following files {} were rejected:\n \n""".format(len(rejected))

    if len(rejected) == 0:
        summary += '    None\n'
    for i in rejected:
        summary += ('    '+i+'\n')

    summary += \
    """    -----------------------------------------------------------------------------
    The following files {} had no status:\n \n""".format(len(no_status))

    if len(accepted) == 0:
        summary += '    None\n'
    for i in no_status:
        summary += ('    '+i+'\n')

    summary += \
    """    =============================================================================
    """

    # Save dataframe and error summary to output.
    df.to_csv('../data/microscopy_status.csv', index=False)

    if os.path.isdir('../data/validation') == False:
        os.mkdir('../data/validation')
    with open('{}{}'.format('../data/validation/', 'microscopy_status_summary.txt'), 'w') as f:
        f.write(summary)

    print(summary)

else:
    print('No dataframes were compiled.')