# -*- coding: utf-8 -*-
# %%
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import skimage.io
import skimage.morphology
import skimage.segmentation
import skimage.filters
import scipy.ndimage
import tqdm
sys.path.insert(0, '../../../')
import mwc.image
import mwc.viz
colors = mwc.viz.personal_style()

# Define teh experimental parameters
DATE = 20180715
TEMP = 37 # in C
CARBON = 'glucose'
OPERATOR = 'O2'
IP_DIST = 0.16

data_dir = '../../../data/images/{}_{}C_{}_{}_bleaching/'.format(DATE, TEMP, CARBON, OPERATOR)
# %% Load the flatfield images and format
noise_files = glob.glob('{}*camera_noise*/*.tif'.format(data_dir))
slide_files = glob.glob('{}/*slide*/*.tif'.format(data_dir))
noise_avgs = {}
field_avgs = {}

# Generate noise images
for i, f in enumerate(noise_files):
    im = skimage.io.imread(f)
    mch = im[0, :, :]
    yfp = im[1, :, :]

    if i == 0:
        noise_avgs['mch'] = mch
        noise_avgs['yfp'] = yfp
    else:
        noise_avgs['mch'] = (noise_avgs['mch'] + mch) / (i + 1)
        noise_avgs['yfp'] = (noise_avgs['yfp'] + yfp) / (i + 1)

# Generate profile images
for i, f in enumerate(slide_files):
    im = skimage.io.imread(f)
    mch = im[0, :, :]
    yfp = im[1, :, :]

    if i == 0:
        field_avgs['mch'] = mch
        field_avgs['yfp'] = yfp
    else:
        field_avgs['mch'] = (field_avgs['mch'] + mch) / (i + 1)
        field_avgs['yfp'] = (field_avgs['yfp'] + yfp) / (i + 1)

# %% Grab the file names of the bleaching frames.
bleaching_files = glob.glob('../../../data/images/{}_{}C_{}_{}_bleaching/photobleaching/*/*.tif'.format(
                DATE, TEMP, CARBON, OPERATOR))

# Instantiate the bleaching data frame.
bleaching_df = pd.DataFrame([])

selem = skimage.morphology.square(3)
cell_id = 0
for i, f in enumerate(bleaching_files):
    im = skimage.io.imread(f)

    # Extract the fluorescence images
    mch_ims = im[:, 1, :, :]

    # Segment the phase image. 
    thresh = skimage.filters.threshold_yen(ims[0, 0, :, :])
    seg = ims[0, 0, :, :] < thresh

    # Prune the objects
    seg = skimage.morphology.remove_small_objects(seg, min_size=100)
    seg = skimage.segmentation.clear_border(seg)
    inv_seg = seg < 1
    seg, num_obj = skimage.measure.label(seg, return_num=True)

    # Process the data
    for j, img in enumerate(tqdm.tqdm(mch_ims)):
        im_flat = mwc.image.generate_flatfield(img, noise_avgs['mch'], field_avgs['mch'])
        mean_bg = np.mean(im_flat[(im_flat * inv_seg) > 0 ])
        props = skimage.measure.regionprops(seg, im_flat)
        for p in props:
            area = p.area * IP_DIST**2
            mean_int = p.mean_intensity
            corrected_intensity = p.area * (mean_int - mean_bg)
            cell_dict = dict(cell_id=cell_id + p.label, corrected_intensity=corrected_intensity,
            area=area, mean_int=mean_int, mean_bg=mean_bg, time=j)
            bleaching_df = bleaching_df.append(cell_dict, ignore_index=True)
    cell_id += num_obj

# %% Process the data
# Compute the mean bleaching curve.
mean_bleach = bleaching_df.groupby(['time']).mean()
mean_norm = mean_bleach['corrected_intensity'] / mean_bleach.iloc[0]['corrected_intensity']
dfs = []

# Normalize the data and compute the residuals.
grouped = bleaching_df.groupby('cell_id')
fig, ax = plt.subplots(1, 1)
for g, d in grouped:
    d = d.copy()
    d.sort_values('time', inplace=True)
    d['norm'] = d['corrected_intensity'] / d.iloc[0]['corrected_intensity']
    d['resid'] = (d['norm'] - mean_norm)**2 * d.iloc[0]['corrected_intensity']

    #Keep only the data greater than zero
    if (d['norm'] <= 1).all() & (d['norm'] >= 0).all()  & (d.iloc[-1]['norm'] < 0.25):
        dfs.append(d)
        _ = ax.plot(d['time'], d['norm'], '-', lw=1, alpha=0.5)
# %%
bins = np.linspace(0, 1, 100)
norm_df = pd.concat(dfs)
norm_df['bin'] = pd.cut(norm_df['norm'], bins)
fig, ax = plt.subplots()
grouped = norm_df.groupby(['bin']).mean()
ax.plot(1 - grouped['norm'], grouped['resid'], 'o')
# norm_df.to_csv('output/{}_{}C_{}_{}_photobleaching.csv'.format(DATE, TEMP, CARBON, OPERATOR))

norm_df

