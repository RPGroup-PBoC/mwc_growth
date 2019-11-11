#%%
import numpy as np 
import pandas as pd
import mwc.bayes

# Load the data 
data = pd.read_csv('../../data/analyzed_foldchange.csv')
data = data[(data['temp']==37) & (data['carbon']=='glucose') & 
            (data['fold_change'] >= 0) & (data['repressors'] > 0)]
data = data.groupby(['date', 'run_number', 'atc_ngml']).mean().reset_index()


# %%
model = mwc.bayes.StanModel('../stan/correction_factor.stan', force_compile=True)

# %%
# Define the data dictionary and sample. 
data_dict = {'N':len(data), 'epRA':-13.9, 'epAI':4.5, 'Nns':4.6E6,
             'foldchange':data['fold_change'], 'R':data['repressors']/2}
fit, samples = model.sample(data_dict)


# %%
