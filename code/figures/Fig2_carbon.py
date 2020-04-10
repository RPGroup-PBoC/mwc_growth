#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob
import skimage.io
import scipy.io
import mwc.viz
import mwc.model
import mwc.process
colors, palette = mwc.viz.plotting_style()

# %%
# Load datasets for data collapse 
old = pd.read_csv('../../data/Garcia2011_Brewster2014_RazoMejia2018_Chure2019.csv')

# Load measurements and restrict
data = pd.read_csv('../../data/analyzed_foldchange.csv', comment='#')
data = mwc.process.condition_filter(data, temp=37)
data['repressors'] *= 1.16

# Compute the predicted bohr parameter
data = data.groupby(['date', 'run_number', 'atc_ngml', 'carbon']).mean()
data = data.groupby(['atc_ngml', 'carbon'])[['fold_change', 'repressors']].agg(('mean', 'sem')).reset_index()

data['bohr'] = -mwc.model.SimpleRepression(R=data['repressors']['mean'], ep_r=-13.9,
                                          ka=139, ki=0.53, ep_ai=4.5, effector_conc=0).bohr_parameter()


# Load the sampling statistics
samps = pd.read_csv('../../data/DNA_binding_energy_samples.csv', comment='#')

# %%

fig, ax  = plt.subplots(2, 2, figsize=(6,4))
ax = ax.ravel()
mwc.viz.despine(ax)
ax[0].set_visible(False)
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlim([1, 300])
ax[1].set_ylim([0.01, 1.1])
ax[3].set_yscale('log')
ax[3].set_xlim([-5, 5])
ax[3].set_ylim([5E-3, 1.5])

# Add labels
ax[1].set_xlabel('repressors per cell', fontsize=8)
ax[1].set_ylabel('fold-change')
ax[2].set_xlabel(r'effective binding energy $\tilde{\epsilon}$ [$k_BT$]')
ax[2].set_ylabel(r'posterior probability')
ax[2].set_yticks([])
ax[3].set_xlabel('free energy [$k_BT$]')
ax[3].set_ylabel('fold-change')

# Plot the theory. 
ep_r = np.array([-14.1, -13.7])
rep_range = np.logspace(0, 3, 300)
ep, r = np.meshgrid(ep_r, rep_range)
fc = mwc.model.SimpleRepression(r, ep, ka=139, ki=0.53, ep_ai=4.5,
                                effector_conc=0).fold_change()

ax[1].fill_between(rep_range, fc[:, 0], fc[:, 1], color=colors['light_grey'], 
                    alpha=0.3, label='theoretical prediction')
carbon_colors = {'glucose':colors['purple'],
                 'glycerol':colors['green'],
                 'acetate':colors['brown']}
for g, d in data.groupby(['carbon']):
    ax[1].errorbar(d['repressors']['mean'], d['fold_change']['mean'], 
                   xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'],
                   markeredgecolor='white', markeredgewidth=0.5,
                   fmt='o', ms=3.75, color=carbon_colors[g], lw=0.5,
                   label=g)
    ax[3].errorbar(d['bohr'], d['fold_change']['mean'], 
                    yerr=d['fold_change']['sem'],
                   markeredgecolor='white', markeredgewidth=0.5,
                   fmt='o', ms=3.75, color=carbon_colors[g], lw=0.5, zorder=1000,
                   label=g)


bohr_range = np.linspace(-10, 10, 200)
master = (1 + np.exp(-bohr_range))**-1
ax[3].plot(old['bohr_parameter'], old['mean'], '.', color='lightgrey', ms=5, 
           label='previous data')
ax[3].plot(bohr_range, master, 'k-', lw=0.75, label='scaling function')

for g, d in samps[samps['temp']==37].groupby(['carbon']):
    bins = np.linspace(-15.1, -13.5, 50)
    ax[2].hist(d[d['parameter']=='epRA']['value'], bins=bins, color=carbon_colors[g],
            edgecolor='white', alpha=0.5, density=True, label=g)


# Add legends and panel labels
ax[1].legend(fontsize=6)
ax[2].legend(fontsize=6)
ax[3].legend(fontsize=6, handlelength=1)
plt.tight_layout()
fig.text(0, 1, '(A)', fontsize=8)
fig.text(0.48, 1, '(B)', fontsize=8)
fig.text(0, 0.48, '(C)', fontsize=8)
fig.text(0.48, 0.48, '(D)', fontsize=8)
plt.savefig('../../figs/Fig2_carbon_sources_plots.pdf', bbox_inches='tight')

# %%
# Load the cell files. 
path = '/Users/gchure/Desktop/dilution_03ngml/xy5'
files = glob.glob(f'{path}/cell/*.mat')

mat = scipy.io.loadmat(files[0])
mat['CellA']

# %%jw
import skimage.measure
import scipy.ndimage
phases = []
yfps = []
mchs = []
for f in files:
    mat = scipy.io.loadmat(f)
    _im = mat['CellA'][0][0]
    mask = _im['mask'][0][0]
    phase = _im['phase'][0][0]
    yfp = _im['fluor1'][0][0]
    mch = _im['fluor2'][0][0]
    props = skimage.measure.regionprops(mask)
    o = -props[0]['orientation'] * (180 / np.pi) + 90
    seg = scipy.ndimage.rotate(mask, o)
    props = skimage.measure.regionprops(seg)
    bbox = props[0]['bbox']
    s = np.s_[bbox[0]:bbox[2], bbox[1]:bbox[3]]
    conts = skimage.measure.find_contours(seg[s], level=0)
    xs = [c[0] for c in conts[0]]
    ys = [c[1] for c in conts[0]]
    phases.append({'x':xs, 'y':ys})
    yfps.append(scipy.ndimage.rotate((yfp * mask), o)[s])
    mchs.append(scipy.ndimage.rotate((mch * mask), o)[s])


#%%
fig, ax = plt.subplots(2, 10, figsize=(2.5, 3), sharex=True)
# Define the order (determined manualy)
order = [0, 1, 6, 8, 3, 4, 5, 9, 2, 7, ]
for i, o in enumerate(order):
    ax[0, i].imshow(yfps[o], cmap='viridis', vmin=360, vmax=500)

    ax[1, i].imshow(mchs[o], cmap='magma', vmin=200, vmax=400, origin='upper')


for i, a in enumerate(ax.ravel()):
    a.set_xlim([0, 15])
    a.set_ylim([0, 60])
    a.set_xticks([])
    a.set_yticks([])
    a.spines['bottom'].set_visible(False)
    a.spines['left'].set_visible(False)
plt.subplots_adjust(hspace=-0.4, wspace=0)
plt.savefig('../../figs/segs.svg', bbox_inches='tight')
# plt.tight_layout()
#  %%


# %%
