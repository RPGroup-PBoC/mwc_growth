import os
import glob
import tqdm
# Our numerical workhorses
import numpy as np
import pandas as pd
import scipy

# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Seaborn, useful for graphics
import seaborn as sns

# Import the project utils
import sys
sys.path.insert(0, '../../../')

import mwc_growth as mwc

# define variables to use over the script
DATE = 20170814
USERNAME = 'gchure'
RUN_NO = 'r1'

# ---------------------------------------------------------------------------
# Shouldn't need to touch anything below here.
# ---------------------------------------------------------------------------

# list the directory with the data
datadir = '../../../data/flow/csv/{0}'.format(DATE)

# Grab all of the files.
files = glob.glob('{0}/*{1}*.csv'.format(datadir, RUN_NO))


# define the parameter alpha for the automatic gating
alpha = 0.40

#  Initialize the DataFrame to save the mean expression levels
cols = ['date', 'username', 'operator', 'temp', 'carbon', 'run_no', 'strain',
        'atc_ngml', 'iptg_um', 'mean_YFP_H']
df = pd.DataFrame([], columns=cols)
# Read the files and compute the mean YFP value
for i, f in tqdm.tqdm(enumerate(files), desc="processing"):
    # Split the file name to extract the information.
    date, operator, temp, carbon, run_no,\
        strain, atc, _, iptg, _ = f.split('/')[-1].split('_')
    temp = int(temp)
    atc = float(atc[:-4])
    iptg = float(iptg[:-2])
    _df = pd.read_csv(f)
    if len(_df) == 0:
        print("No data in file {0}".format(f))
        pass
    else:
        # Gate on an arbitrary percentile of a bivariate lognormal distribution
        # of forward scattering and side scattering.
        gated = mwc.auto_gauss_gate(_df, alpha, x_val='FSC-A',
                                    y_val='SSC-A', log=True)

        # Compute the mean YFP value
        mean_YFP = gated['FITC-H'].mean()
        df_vars = [date, USERNAME, operator,
                   temp, carbon, run_no, strain, atc, iptg, mean_YFP]
        # Read the file.
        df_dict = {v: [df_vars[i]] for i, v in enumerate(cols)}
        _df = pd.DataFrame(df_dict)
        df = df.append(_df, ignore_index=True)

# Save the dataframe to disk.
target = 'output/{0}_{1}_{2}_{3}_{4}_{5}_IPTG_titration_MACSQuant.csv'.format( DATE, RUN_NO, USERNAME, operator, temp, carbon)
if os.path.exists(target) is True:
    os.remove(target)
with open('comments.txt', 'r') as f:
    comments = f.readlines()
with open(target, 'a') as f:
    for line in comments:
        f.write(line)
    df.to_csv(f, index=False)
