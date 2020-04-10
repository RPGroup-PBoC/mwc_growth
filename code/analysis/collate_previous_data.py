#%%
import numpy as np
import pandas as pd
import mwc.model
import mwc.viz
colors, palette = mwc.viz.plotting_style()
constants = mwc.model.load_constants()

# Load the mutant sampling information
mutants_epr = pd.read_csv('../../data/Chure2019_DNA_binding_energy_summary.csv')
mutants_allo = pd.read_csv('../../data/Chure2019_KaKi_epAI_summary.csv')

# Load and clean the datasets
induction = pd.read_csv('../../data/RazoMejia_2018.csv', comment='#')
induction = induction[induction['repressors'] > 0]
induction['repressors'] *= 2
induction = induction.groupby(['IPTG_uM', 'repressors', 
                                'operator'])['fold_change_A'].agg(
                                    ('mean', 'sem')).reset_index()
induction.rename(columns={'fold_change_A':'fold_change', 'IPTG_uM':'IPTGuM'},
                 inplace=True)
mutants = pd.read_csv('../../data/Chure2019_summarized_data.csv', comment='#')
garcia_brewster = pd.read_csv('../../data/Garcia2011_Brewster2014.csv', comment='#')
garcia_brewster['IPTGuM'] = 0
garcia_brewster.rename(columns={'repressor':'repressors', 'fold_change':'mean'}, inplace=True)

# Compute the bohr parameters
induction['bohr_parameter'] = -mwc.model.SimpleRepression(induction['repressors'].values, 
                              ep_r=[constants[o] for o in induction['operator'].values],
                              ka=constants['ka'], ki=constants['ki'], 
                              ep_ai=constants['ep_ai'], 
                              effector_conc=induction['IPTGuM'].values
                              ).bohr_parameter()
garcia_brewster['bohr_parameter'] = -mwc.model.SimpleRepression(garcia_brewster['repressors'].values, 
                              ep_r=[constants[o] for o in garcia_brewster['operator'].values],
                              ka=constants['ka'], ki=constants['ki'], 
                              ep_ai=constants['ep_ai'], 
                              effector_conc=0).bohr_parameter()
garcia_brewster['mutant'] = 'wt'
# Mutants is more complicated.
DNA  = {'Y20I':-9.9, 'Q21M':-15.4, 'Q21A':-11}
IND = {'F164T': {'ka':165, 'ki':3, 'ep_ai':1},
       'Q294K': {'ka':1000, 'ki':310, 'ep_ai':-3.11},
       'Q294R': {'ka':9, 'ki':8, 'ep_ai':-2.6},
       'Q294V': {'ka':650, 'ki':8, 'ep_ai':0.1}}
dfs = []
for g, d in mutants.groupby(['class', 'mutant', 'operator', 'repressors']):
    if g[0] == 'DNA':
        ep_r = mutants_epr[(mutants_epr['parameter']=='ep_RA') & 
                           (mutants_epr['mutant']==g[1]) & 
                           (mutants_epr['repressors']==260)]['median'].values[0]
        ka = constants['ka']
        ki = constants['ki']
        ep_ai = constants['ep_ai']
    elif g[0] == 'IND':
        allo = mutants_allo[(mutants_allo['mutant']==g[1]) & 
                            (mutants_allo['operator']=='O2')]
        ep_r = np.array([constants[o] for o in d['operator'].values])
        ka = allo[allo['parameter']=='Ka']['median'].values[0]
        ki = allo[allo['parameter']=='Ki']['median'].values[0]
        ep_ai = allo[allo['parameter']=='ep_AI']['median'].values[0]
    elif g[0] == 'DBL':

        _DNA, _IND = g[1].split('-')
        ep_r = mutants_epr[(mutants_epr['parameter']=='ep_RA') & 
                           (mutants_epr['mutant']==_DNA) & 
                           (mutants_epr['repressors']==260)]['median'].values[0]
        allo = mutants_allo[(mutants_allo['mutant']==_IND) & 
                            (mutants_allo['operator']=='O2')]
        ka = allo[allo['parameter']=='Ka']['median'].values[0]
        ki = allo[allo['parameter']=='Ki']['median'].values[0]
        ep_ai = allo[allo['parameter']=='ep_AI']['median'].values[0]
    elif g[0] == 'WT':
        ep_r = np.array([constants[o] for o in d['operator'].values])
        ka = constants['ka']
        ki = constants['ki']
        ep_ai = constants['ep_ai']
    d = d.copy()
    d['bohr_parameter'] = -mwc.model.SimpleRepression(R=d['repressors'].values,
                           ep_r=ep_r, ka=ka, ki=ki, ep_ai=ep_ai, 
                           effector_conc=d['IPTGuM'].values).bohr_parameter()
    d['author'] = 'chure'
    d = d[['author', 'operator', 'mean', 'sem', 'mutant', 'IPTGuM', 
          'repressors', 'bohr_parameter']]
    dfs.append(d) 
mutants = pd.concat(dfs)

# concatenate and save
induction['author'] = 'razo-mejia'
induction['mutant'] = 'wt'
wt_rep = pd.concat([induction, garcia_brewster, mutants], sort=False)
wt_rep.drop(columns=['energy'], inplace=True)
wt_rep.to_csv('../../data/Garcia2011_Brewster2014_RazoMejia2018_Chure2019.csv', index=False)




# %%
