#%%
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import mwc.bayes
import mwc.stats
import mwc.model
import matplotlib.transforms
colors, _ = mwc.viz.personal_style()

# %%
# Load the various data sets. 
data = pd.read_csv('../../data/analyzed_foldchange.csv')
data = data[(data['carbon']=='glucose') & (data['temp']==37) & 
            (data['strain']=='dilution') & (data['repressors'] > 0) & 
            (data['fold_change'] >= 0)]

# Summarize the data
large_only = data[data['size']=='large']
data = data.groupby(['date', 'run_number', 'atc_ngml']).mean().reset_index()
large_only = large_only.groupby(['date', 'run_number', 'atc_ngml']).mean().reset_index()


# Load the garcia and brewster O2 data
old_gods = pd.read_csv('../../data/Garcia2011_Brewster2014.csv', comment='#')
old_gods = old_gods[old_gods['operator']=='O2']
old_gods.rename(columns={'repressor':'repressors'}, inplace=True)

garcia = old_gods[old_gods['author']=='garcia']
brewster = old_gods[old_gods['author']=='brewster']

# Load the Razo-Mejia O2 data. 
ind_data = pd.read_csv('../../data/RazoMejia_2018.csv', comment='#')
ind_data.rename(columns={'fold_change_A':'fold_change'}, inplace=True)
ind_data['repressors'] *= 2
ind_data = ind_data[(ind_data['operator']=='O2') & (ind_data['repressors'] > 0) &
                    (ind_data['IPTG_uM']==0)]
# %%
# Load the stan model
model = mwc.bayes.StanModel('../stan/DNA_binding_energy.stan')

# %%
# Perform the inference. 
summ_dfs = []
data_dict = {'no_correction':data, 'correction':data, 'large_only':large_only, 
            'all_divided':data, 'garcia':garcia,
            'brewster':brewster, 'razo-mejia':ind_data}
for k, v in data_dict.items():
    # Define the data dictionary. 
    data_dict = {'N':len(v), 'foldchange':v['fold_change'],
                'Nns':4.6E6}
    if k == 'no_correction':
        data_dict['repressors'] = v['raw_repressors']

    elif k == 'all_divided':
        data_dict['repressors'] = v['raw_repressors'].values * 2

    else:
        data_dict['repressors'] = v['repressors']
    fit, samples = model.sample(data_dict)
    params = model.summarize_parameters()
    params['source'] = k 
    summ_dfs.append(params)
stats = pd.concat(summ_dfs)

# Save the sampling summary to disk. 
stats.to_csv('../../data/R_correction_DNA_binding_energy_summary.csv',
            index=False)

# # %%
# # Set up a figure to plot the binding energy values
# fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=100)

# # Adjust the axes. 
# ax.set_ylim([-16, -12])
# ax.set_xlabel('data set')
# ax.set_ylabel('DNA binding energy [$k_BT$]')
# mwc.viz.titlebox(ax, r'operator O2, $\Delta\varepsilon_{AI} \rightarrow \infty$', color=colors['black'], boxsize="15%")

# position = {'present work':0, 'Garcia & Phillips 2011':1, 
#            'Brewster et al. 2014':2, 'Razo-Mejia et al. 2018': 3, 'pooled':4}
# edge_colors = {'present work':colors['dark_purple'], 
#                'Garcia & Phillips 2011':colors['dark_orange'], 
#                 'Brewster et al. 2014':colors['dark_blue'], 
#                 'Razo-Mejia et al. 2018': colors['dark_red'],'pooled':colors['black']}

# face_colors = {'present work':colors['light_purple'], 
#                'Garcia & Phillips 2011':colors['light_orange'], 
#                 'Brewster et al. 2014':colors['light_blue'], 
#                 'Razo-Mejia et al. 2018': colors['light_red'], 'pooled':colors['black']}


# for g, d in stats.groupby('source'):
#     median, low, high = d[d['parameter']=='epRA'][['median', 
#                                                 'hpd_min', 'hpd_max']].values[0]
#     ax.vlines(position[g], low, high, lw=1, color=edge_colors[g])
#     ax.plot(position[g], median, '.', ms=8, markeredgecolor=edge_colors[g],
#             markerfacecolor=face_colors[g], markeredgewidth=0.75, label=g)

# ax.set_xticks([v for _, v in position.items()])
# ax.set_xticklabels([k for k, _ in position.items()])
# plt.savefig('../../figs/FigSX_DNA_binding_energy_comparison.pdf', 
#             bbox_inches='tight', dpi=300, facecolor='white')
# # %%
# # Generate the pairwise plot.
# fig, ax = plt.subplots(4, 4, sharex=True, sharey=True, figsize=(6.5, 6.5), dpi=100)
# for a in ax.ravel():
#     a.set_xscale('log')
#     a.set_yscale('log')
#     a.set_xlim([1, 1E4])
#     a.set_ylim([1E-3, 1.2])
#     a.set_yticks([1, 1E-1, 1E-2])
#     a.set_xticks([10, 100, 1000])

# # Change the diagonals    
# for i in range(4):
#     ax[i, i].set_facecolor('white')
#     ax[i, i].grid(color=colors['grey'])
#     ax[-1, i].set_xlabel('rep. / cell', fontsize=8)
#     ax[i, 0].set_ylabel('fold-change', fontsize=8)

# for k, v in position.items():
#     if k != 'pooled':
#         mwc.viz.ylabelbox(ax[v, 0], k, color=edge_colors[k], boxsize="20%", fontsize=6)
#         mwc.viz.titlebox(ax[0, v], k, color=edge_colors[k], boxsize="20%", fontsize=6)
#         dx = -21 / fig.dpi
#         dy = 0
#         offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
#         if v !=  0:
#             for label in ax[v, 0].yaxis.get_majorticklabels():
#                 label.set_transform(label.get_transform() + offset)

# # Plot the points in each. 
# for d, s in zip(desc[:-1], src[:-1]):
#     for i in range(4):
#         pos = position[d]
#         if d == 'present work':
#             _s = s.groupby(['atc_ngml']).agg(('mean', 'sem')).reset_index()     
#             ax[i, pos].errorbar(_s['repressors']['mean'], _s['fold_change']['mean'],
#                             xerr=_s['repressors']['sem'], 
#                             yerr=_s['fold_change']['sem'], color=edge_colors[d],
#                             fmt='.', ms=8, markeredgecolor=edge_colors[d], 
#                             markerfacecolor=face_colors[d], markeredgewidth=0.75, 
#                             lw=0.75, capsize=1)
#         elif d == 'Razo-Mejia et al. 2018':
#             _s = s.groupby(['repressors']).agg(('mean', 'sem')).reset_index()
#             ax[i, pos].errorbar(_s['repressors'], _s['fold_change']['mean'],
#                             yerr=_s['fold_change']['sem'], color=edge_colors[d],
#                             fmt='.', ms=8, markeredgecolor=edge_colors[d], 
#                             markerfacecolor=face_colors[d], markeredgewidth=0.75, 
#                             lw=0.75, capsize=1)

#         else:
#             ax[i, pos].plot(s['repressors'], s['fold_change'], '.', ms=8, 
#                            markeredgecolor=edge_colors[d], 
#                            markerfacecolor=face_colors[d], markeredgewidth=0.75)


# # Plot the predictions/fits
# rep_range = np.logspace(0, 4, 200)
# for g, d in stats.groupby(['source']):
#     low, high= d[d['parameter']=='epRA'][['hpd_min', 'hpd_max']].values[0]
#     low_fit = mwc.model.SimpleRepression(R=rep_range, ep_r=low, ep_ai=4.5, ka=139,
#                                     ki=0.53, n_sites=2, n_ns=4.6E6, effector_conc=0).fold_change()
#     high_fit = mwc.model.SimpleRepression(R=rep_range, ep_r=high, ep_ai=4.5, ka=139,
#                                     ki=0.53, n_sites=2, n_ns=4.6E6, effector_conc=0).fold_change()
#     if g != 'pooled':
#         for i in range(4):
#             ax[position[g], i].fill_between(rep_range, low_fit, high_fit, 
#                         color=face_colors[g], alpha=0.5) 
#     # else:
#     #     for i in range(4):
#     #         for j in range(4):
#     #             ax[i,j].fill_between(rep_range, low_fit, high_fit, 
#     #                     color=face_colors[g], alpha=0.15) 


# plt.savefig('../../figs/FigSX_pairwise_DNA_binding_energy_comparison.pdf', 
#             bbox_inches='tight', dpi=100, facecolor='white')
# plt.subplots_adjust(wspace=0.05, hspace=0.05)


# # %%
# data = pd.read_csv('../../data/analyzed_foldchange.csv')
# data = data[(data['date'] != 20190307) & (data['fold_change'] >= 0) & 
#             (data['repressors'] > 0)]

# # Load the cell size distributions. 
# sizes = pd.read_csv('../../data/analyzed_fluctuations.csv')

# data['repressors'] /= 2
# carbon = 'glucose'
# temp = 42 
# dfs = []
# for g, d in data.groupby(['date', 'run_number', 'carbon', 'temp', 'strain', 'atc_ngml']):
#     d = d.copy()
#     birth_length = sizes[(sizes['date']==int(g[0])) & (sizes['run_no']==g[1]) & 
#                          (sizes['carbon']==g[2]) & (sizes['temp']==g[3]) 
#                          ][['length_1_birth', 'length_2_birth']].values.flatten() / 0.065
#     thresh_ind = np.where(np.linspace(0, 1, len(birth_length)) >=0.90)[0][0]
#     thresh = np.sort(birth_length)[thresh_ind]
#     d.loc[d['volume_birth'] < thresh, 'size'] = 'small'
#     d.loc[d['volume_birth'] >= thresh, 'size'] = 'large'
#     d.loc[d['size']=='small', 'repressors'] *= 2
#     dfs.append(d)
# data = pd.concat(dfs)

# # fig, ax = plt.subplots(dpi=100)
# # _ = ax.hist(data['length_um'], bins=bins, alpha=0.5, label='data', density=True)
# # _ = ax.hist(birth_length, bins=bins, alpha=0.5, label='birth', density=True)
# # _ = ax.vlines(thresh, 0, ax.get_ylim()[1], color=colors['red'], lw=2, label=thresh)
# # ax.legend()

# data = data[(data['strain']=='dilution') & 
#             (data['carbon']==carbon) & (data['temp']==temp)]
# grouped = data.groupby(['date', 'run_number', 'atc_ngml']).mean().reset_index()
# grouped = grouped.groupby(['atc_ngml']).agg(('mean', 'sem')).reset_index()

# fig, ax = plt.subplots(1, 1, dpi=100)
# ax.set_xscale('log')
# rep_range = np.logspace(0, 3, 200)
# theo = mwc.model.SimpleRepression(R=rep_range, ep_r=-13.9, effector_conc=0, ka=139,
#                                 ki=0.53, ep_ai=4.5).fold_change()
# ax.plot(rep_range, theo)
# ax.errorbar(grouped['repressors']['mean'], grouped['fold_change']['mean'], 
#             yerr=grouped['fold_change']['sem'], xerr=grouped['repressors']['sem'], fmt='.')

# ax.legend()
# ax.set_yscale('log')
# #%%
# # Redo the inference for the glucose case. 
# # %%
# # Perform the inference. 

# summ_dfs = []
# desc = ['present work', 'Garcia & Phillips 2011', 'Brewster et al. 2014',
#         'Razo-Mejia et al. 2018', 'pooled']
# _data = data.groupby(['date', 'run_number', 'atc_ngml']).mean()
# src = [_data, garcia, brewster, ind_data]
# src.append(pd.concat(src))
# for d, s  in zip(desc, src):
#     # Define the data dictionary. 
#     data_dict = {'N':len(s), 'foldchange':s['fold_change'],
#                 'repressors':s['repressors'], 'Nns':4.6E6}
#     fit, samples = model.sample(data_dict)
#     params = model.summarize_parameters()
#     params['source'] = d
#     summ_dfs.append(params)

# refit_stats = pd.concat(summ_dfs)


# # %%
# # Set up a figure to plot the binding energy values
# fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=100)

# # Adjust the axes. 
# ax.set_ylim([-16, -12])
# ax.set_xlabel('data set')
# ax.set_ylabel('DNA binding energy [$k_BT$]')
# mwc.viz.titlebox(ax, r'operator O2, $\Delta\varepsilon_{AI} \rightarrow \infty$', color=colors['black'], boxsize="15%")

# position = {'present work':0, 'Garcia & Phillips 2011':1, 
#            'Brewster et al. 2014':2, 'Razo-Mejia et al. 2018': 3, 'pooled':4}
# edge_colors = {'present work':colors['dark_purple'], 
#                'Garcia & Phillips 2011':colors['dark_orange'], 
#                 'Brewster et al. 2014':colors['dark_blue'], 
#                 'Razo-Mejia et al. 2018': colors['dark_red'],'pooled':colors['black']}

# face_colors = {'present work':colors['light_purple'], 
#                'Garcia & Phillips 2011':colors['light_orange'], 
#                 'Brewster et al. 2014':colors['light_blue'], 
#                 'Razo-Mejia et al. 2018': colors['light_red'], 'pooled':colors['black']}


# for g, d in refit_stats.groupby('source'):
#     median, low, high = d[d['parameter']=='epRA'][['median', 
#                                                 'hpd_min', 'hpd_max']].values[0]
#     ax.vlines(position[g], low, high, lw=1, color=edge_colors[g])
#     ax.plot(position[g], median, '.', ms=8, markeredgecolor=edge_colors[g],
#             markerfacecolor=face_colors[g], markeredgewidth=0.75, label=g)

# ax.set_xticks([v for _, v in position.items()])
# ax.set_xticklabels([k for k, _ in position.items()])
# plt.savefig('../../figs/FigSX_refit_DNA_binding_energy_comparison.pdf', 
#             bbox_inches='tight', dpi=300, facecolor='white')
# #

# # %%


# %%
