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
import skimage.measure
import scipy.ndimage
import scipy.signal
colors, palette = mwc.viz.plotting_style()

# #############################################################################
# Generate plots of fold-change and binding energy distributions.
# ############################################################################ 
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

fig, ax  = plt.subplots(2, 2, figsize=(6,5))
ax = ax.ravel()
mwc.viz.despine(ax)
ax[0].set_visible(False)
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlim([1, 300])
ax[1].set_ylim([0.005, 1.1])
ax[3].set_yscale('log')
ax[3].set_xlim([-5, 5])
ax[3].set_ylim([5E-3, 1.5])

# Add labels
ax[1].set_xlabel('repressors per cell', fontsize=8)
ax[1].set_ylabel('fold-change')
ax[2].set_xlabel(r'effective binding energy $\epsilon$ [$k_BT$]')
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
                   markeredgecolor=colors['grey'], markeredgewidth=0.5,
                   fmt='o', ms=3.75, color=carbon_colors[g], lw=0.5,
                   label=g)
    ax[3].errorbar(d['bohr'], d['fold_change']['mean'], 
                    yerr=d['fold_change']['sem'],
                   markeredgecolor=colors['grey'], markeredgewidth=0.5,
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
            edgecolor=colors['grey'], alpha=0.5, density=True, label=g, lw=0.25)


# Add legends and panel labels
ax[1].legend(fontsize=6)
ax[2].legend(fontsize=6)
ax[3].legend(fontsize=6, handlelength=1)
plt.tight_layout()
plt.savefig('../../figs/Fig2_carbon_sources_plots.pdf', bbox_inches='tight')

#%%
# #############################################################################
# Generate the image plot
# ############################################################################

# %%
# Load the cell files. 
files = glob.glob(f'../../data/example_cells/*.mat')

# Instantiate empty lists to store the cell images and contours
contours = []
yfps = []
mchs = []

# Iterate through each cell file
for f in files:

    # Load the mat and separate the images
    mat = scipy.io.loadmat(f)
    _im = mat['CellA'][0][0]
    mask = _im['mask'][0][0]
    yfp = _im['fluor1'][0][0]
    mch = _im['fluor2'][0][0]

    # From the segmentation mask, compute the properties and the rotation 
    props = skimage.measure.regionprops(mask)
    o = -props[0]['orientation'] * (180 / np.pi)

    # Rotate the segmentation mask and find the bounding box
    seg = scipy.ndimage.rotate(mask, o)
    props = skimage.measure.regionprops(seg)
    bbox = props[0]['bbox']

    # Pad the bounding box by 10 px
    s = np.s_[bbox[0]-5:bbox[2]+5, bbox[1]-5:bbox[3]+5]

    # From the segmentation mask, compute the segmentation contours
    conts = skimage.measure.find_contours(seg[s], level=0)
    xs = [c[0] for c in conts[0]]
    ys = [c[1] for c in conts[0]]
    xs = scipy.signal.savgol_filter(xs, window_length=7, polyorder=2)
    ys = scipy.signal.savgol_filter(ys, window_length=7, polyorder=2)
    # Add the contours and fluorescence images to the storage list
    contours.append({'x':xs, 'y':ys})
    yfps.append(scipy.ndimage.rotate((yfp * mask), o)[s])
    mchs.append(scipy.ndimage.rotate((mch * mask), o)[s])

#%%
fig, ax = plt.subplots(10, 2, figsize=(3, 2.5), sharey=True)
# Define the order (determined manualy)
order = [0, 1, 6, 8, 3, 4, 5, 9, 2, 7, ]
for i, o in enumerate(order):
    ax[i, 0].imshow(yfps[o], cmap='viridis', vmin=365, vmax=500)
    ax[i, 1].imshow(mchs[o], cmap='magma', vmin=210, vmax=380, origin='upper')
    for j in range(2):
        ax[i, j].plot(contours[o]['y'], contours[o]['x'], 'w-', lw=0.5)

# Format the axes so the YFP and mCherry images are completely reflected
for i in range(10):
    ax[i, 0].set_xlim([65, 3])
    ax[i, 0].set_facecolor('#440154')
    ax[i, 1].set_xlim([3, 65])
    ax[i, 1].set_facecolor('#000003')

for i, a in enumerate(ax.ravel()):
    a.set_xticks([])
    a.set_yticks([])
    a.spines['bottom'].set_visible(False)
    a.spines['left'].set_visible(False)
plt.subplots_adjust(hspace=-0.1, wspace=-0.7)
plt.savefig('../../figs/Fig2_segmentation.svg', bbox_inches='tight')
# %%
