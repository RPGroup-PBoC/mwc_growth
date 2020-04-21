#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mwc.viz
import mwc.model
colors, palette = mwc.viz.plotting_style()

# Load the data sets. 
data = pd.read_csv('../../data/Garcia2011_Brewster2014_RazoMejia2018_Chure2019.csv', comment='#')

#%%
fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5))
mwc.viz.despine(ax)
ax.set_xlabel('free energy [$k_BT$]')
ax.set_ylabel('fold-change')
ax.set_yscale('log')
ax.set_xlim([-8, 9])
ax.set_ylim([5E-4, 1.8])

# Define the bohr parameter
bohr_range = np.linspace(-10, 10, 300)
collapse = (1 + np.exp(-bohr_range))**-1

op_glyphs = {'O1':'o', 'O2':'^', 'O3': 'D', 'Oid':'v'}
rep_colors = {22:colors['purple'], 60:colors['red'], 124:colors['green'],
              260:colors['blue'], 1220:colors['brown'], 1740:colors['orange']}

for g, d in data.groupby(['repressors', 'IPTGuM', 
                          'operator', 'author',
                          'mutant']):
    if g[0] not in list(rep_colors.keys()):
        _color = colors['black']
    else:
        _color = rep_colors[g[0]]
    ax.plot(d['bohr_parameter'], d['mean'], marker=op_glyphs[g[2]],
            color=_color, markeredgewidth=0.5, markeredgecolor=colors['grey'],
            linestyle='none', ms=3.5, alpha=0.75)
ax.plot(bohr_range, collapse, 'k-', lw=0.75, label='scaling function')
plt.savefig('../../figs/fig1_collapse.svg', bbox_inches='tight')
# %%

