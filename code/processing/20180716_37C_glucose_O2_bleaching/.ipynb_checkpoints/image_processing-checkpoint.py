# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import skimage.io
import skimage.morphology
import skimage.segmentation
import skimage.filters
import scipy.ndimage
import glob
import matplotlib.pyplot as plt
import tqdm
sys.path.insert(0, '../../../')
import mwc.viz
import mwc.image
import mwc.stats
import holoviews as hv
import bokeh.io
colors = mwc.viz.personal_style()
hv.extension('matplotlib')

# %%
# Define the experimental parameters
DATE = 20180716
CARBON = 'glucose'
TEMP = 37
OPERATOR = 'O2'
IP_DIST = 0.07
num_pos = 8 
data_dir = '../../../data/images/{}_{}C_{}_{}_bleaching'.format(DATE, TEMP, CARBON, OPERATOR)

# Load the slide images 
slide_ims = skimage.io.ImageCollection(glob.glob('{}/*slide*.TIF'.format(data_dir)))
noise_ims = skimage.io.ImageCollection(glob.glob('{}/*noise*.TIF'.format(data_dir)))
field_avg = np.mean(slide_ims, axis=0)
noise_avg = np.mean(noise_ims, axis=0)

# %% Find all bleaching positions. 
bleaching_df = pd.DataFrame([])
cell_id = 0
for i in tqdm.tqdm(range(num_pos)):
    files = np.sort(glob.glob('{}/photobleaching/*xy{:02d}_t*.tif'.format(data_dir, i)))

    # Load and segment
    im_init = skimage.io.imread(files[0])
    thresh = skimage.filters.threshold_li(im_init)
    seg = im_init > thresh
    inv_seg = im_init < thresh
    seg = skimage.morphology.remove_small_objects(seg, min_size=100)
    lab, num_obj = skimage.measure.label(seg, return_num=True)

    # Loop though each image, flatten, and extract info
    for j, f in enumerate(files):
        im = skimage.io.imread(f)
        flat = mwc.image.generate_flatfield(im, noise_avg, field_avg)
        mean_bg = np.mean(flat[(inv_seg * flat) > 0])
        props = skimage.measure.regionprops(lab, flat) 
        for p in props:
            area = p.area * IP_DIST**2
            mean_int = p.mean_intensity
            corrected = (mean_int - mean_bg) * p.area
            cell_dict = dict(cell_id=int(p.label + cell_id), mean_intensity=mean_int, area=area,
                              time=j, corrected_intensity=corrected)
            bleaching_df = bleaching_df.append(cell_dict, ignore_index=True)
    cell_id += num_obj
# %%
# Compute the mean bleaching trace.
mean_bleach = bleaching_df.groupby('time').median()

# Group by cell ID and compute normalized value and  the mean bleaching trend
grouped = bleaching_df.groupby('cell_id')
dfs = []
for g, d in grouped:
    if len(d) == len(mean_bleach):
        d = d.copy()
        d.sort_values('time', inplace=True)
        d['norm'] = (d['corrected_intensity'].values - d['corrected_intensity'].values.min()) / (d['corrected_intensity'].values.max() - d['corrected_intensity'].values.min())
        mean_norm = (mean_bleach['corrected_intensity'].values  - mean_bleach['corrected_intensity'].min())/ (mean_bleach['corrected_intensity'].max() - mean_bleach['corrected_intensity'].min())
        d['resid'] = (d['norm'].values - mean_norm)**2 * d.iloc[0]['corrected_intensity']
        dfs.append(d)
norm_df = pd.concat(dfs)
norm_df.to_csv('output/{}_{}C_{}_{}_bleaching.csv'.format(DATE, TEMP, CARBON, OPERATOR))

# %%
bins = np.linspace(0, 1, 5)
norm_df['bins'] = pd.cut(norm_df['norm'], bins)
grouped = pd.DataFrame(norm_df.groupby('bins').median()).reset_index()
grouped.dropna(inplace=True)
alpha = 2 * 6 * np.trapz(grouped['resid'], grouped['norm'])
hv.Curve(grouped, ['norm'], ['resid'])
norm_df[norm_df['time']==0]['corrected_intensity'].values / alpha



bs = mwc.stats.fast_bootstrap(norm_df, 5, iter=1E4)
hv.Distribution(bs)

# traces = hv.Curve(norm_df, ['time', 'cell_id'], ['norm']).groupby('cell_id').options(alpha=0.05, color='slategray').overlay()
# mean = hv.Curve((np.arange(0, len(mean_norm), 1), mean_norm)).options(color='firebrick')
# traces * mean