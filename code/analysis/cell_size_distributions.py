#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mwc.stats
import mwc.io
import mwc.process
from  tqdm import tqdm
import glob

# Parse through all approved experiments
exps = glob.glob('../../data/preprocessed/20*/growth/xy*/clist.mat')

dfs = []
for i, exp in enumerate(tqdm(exps)):
   df = mwc.process.clist_to_dataframe(exp) 
   date, run_number, temp, carbon, _, _ = exp.split('preprocessed/')[1].split('/')[0].split('_')
   temp = int(temp[:-1])
   run_number = int(run_number[1:])
   df['date'] = date
   df['temp'] = temp
   df['run_number'] = run_number
   df['carbon'] = carbon

   # Remove cells with segmentation errors.
   df['valid'] = df['error_frame'].isnull()
   df = df[df['valid']==True]

   # Remove cells that were born on the first frame or died on the last frame. 
   df = df[(df['cell_birth_time'] != 1) & 
          (df['cell_birth_time'] != df['cell_birth_time'].max())]
   dfs.append(df)
data = pd.concat(dfs)

# %%
data.to_csv('../../data/cell_sizes.csv')

# %%
