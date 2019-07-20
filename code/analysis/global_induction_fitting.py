#%%
import numpy as np
import pandas as pd
import mwc.bayes
import mwc.stats
import mwc.viz
import bokeh.io
import bokeh.plotting
bokeh.io.output_notebook()
colors, color_list = mwc.viz.bokeh_theme()

# Load the Garcia and Induction data. 
garcia = pd.read_csv('../../data/GarciaPhillips_2011.csv', comment='#')
razo = pd.read_csv('../../data/RazoMejia_2018.csv', comment='#')

# Get the mapping of repressor count errors to RBS
errs = {'RBS1':80, 'HG104':2, 'RBS1147':10, 'RBS1027':20, 'RBS446':15, 'RBS1L':170}

# Add the errors to the razo data. 
for r, err in errs.items():
    razo.loc[razo['rbs']==r, 'repressors_err'] = err

# Remove the unknown (to me) RBS from the garcia data. 
garcia = garcia[(garcia['rbs']=='RBS1') | (garcia['rbs']=='HG104') | 
               (garcia['rbs']=='RBS1147') | (garcia['rbs']=='RBS1027') |
               (garcia['rbs']=='RBS446') | (garcia['rbs']=='1I')].copy()

# Drop Oid
garcia = garcia[garcia['operator'] != 'Oid'].copy()

# Rename '1I'to the new name
garcia.loc[garcia['rbs']=='1I', 'rbs'] = 'RBS1L'

# Make sure they all have the same column names. 
garcia.rename(columns={'delta_repressors':'repressors_err'}, inplace=True) 
razo.rename(columns={'fold_change_A':'fold_change', 'IPTG_uM':'IPTGuM'}, inplace=True)

# Add author information and missing fields 
garcia['reference'] = 'garcia'
garcia['IPTGuM'] = 0
razo['reference'] = 'razo-mejia'

# Reduce and concatenate
dfs = [garcia[['rbs', 'operator', 'fold_change', 'repressors', 
               'repressors_err', 'IPTGuM', 'reference']],
      razo[['rbs', 'operator', 'fold_change', 'repressors', 
            'repressors_err', 'IPTGuM', 'reference']]]
data = pd.concat(dfs)

# Get rid of any  measuremetns with zoer repressors per cell
data = data[data['repressors'] > 0]

# Convert tetramer counts to dimers. 
data['repressors'] *= 2
data['repressors_err'] *= 2

# Add unique identifiers for the repressor IDs and oeprator IDs
data['rep_idx'] = data.groupby('rbs').ngroup() + 1
data['op_idx'] = data.groupby('operator').ngroup() + 1
#%%
# Load the stan model. 
model = mwc.bayes.StanModel('../stan/global_glucose_param_estimation.stan',
                           force_compile=True)


#%%
# Assemble the data dictionary. 
data_dict = {'N':len(data), 'NR': np.max(data['rep_idx'].unique()),
             'Nep':np.max(data['op_idx'].unique()), 
             'R_idx':data['rep_idx'], 'ep_idx':data['op_idx'],
             'R_mu':data['repressors'].unique(), 
             'R_sig':data['repressors_err'].unique(),
             'ep_ai':4.5, 'n_sites':2, 'n_ns':4.6E6, 'fc_obs':data['fold_change'],
             'c':data['IPTGuM']}

# Sample!
fit, sample = model.sample(data_dict, iter=5000, control=dict(adapt_delta=0.99))

#%%
p = bokeh.plotting.figure(x_axis_type='log')
p.circle(data['IPTGuM'], data['fold_change'])
bokeh.io.show(p)

#%%
