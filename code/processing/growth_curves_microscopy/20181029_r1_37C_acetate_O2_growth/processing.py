# -*- coding: utf-* -*-
import sys
import os
import numpy as np
import pandas as pd
import glob
import skimage.io
import skimage.measure
import skimage.morphology
import skimage.segmentation
import skimage.filters
import scipy.ndimage
import matplotlib.pyplot as plt
sys.path.insert(0, '../../../../')
import mwc.viz
import mwc.image
colors = mwc.viz.personal_style()

# Define the constants
DATE = 20181029
RUN_NO = 1
CARBON = 'acetate'
TEMP = 37  # in C
OPERATOR = 'O2'
INTERVAL = 10  # in min
IP_DIST = 0.065  # in µm per pixel

# Make the output directory if necessary. 
if os.path.isdir('output') == False:
    os.mkdir('output')

# Determine the max number of poisitions
files = np.sort(glob.glob(f'/Volumes/GDC_DATA_2/{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{OPERATOR}_growth/*c2.tif'))
NUM_POS = int(files[-1].split('xy')[-1].split('c')[0])

# Determine the max number of positions
df = pd.DataFrame([], columns=['date', 'carbon', 'temp_C', 'operator', 'colony_idx', 'area', 'time_min'])
                    
col_no = 1
for i in range(NUM_POS): 
    files = glob.glob(
        '/Volumes/GDC_DATA_2/{}_r{}_{}C_{}_{}_growth/*xy{:02d}c2*.tif'.format(DATE, RUN_NO,
                                                                TEMP, CARBON, 
                                                                OPERATOR, i+1))
    files = np.sort(files)
    ims = skimage.io.ImageCollection(files)
    files

    # Correct for out of register images.
    shifted = mwc.image.correct_drift(ims, verbose=True)
    selem = skimage.morphology.square(5)
    im_filt = [scipy.ndimage.median_filter(i, footprint=selem) for i in shifted]

    # Set the global threshold for the image.
    im_thresh = [im >= skimage.filters.threshold_otsu(im) for im in im_filt]

    # %% identify the bounding boxes of all terminal microcolonies 
    cleared = skimage.segmentation.clear_border(im_thresh[-1], buffer_size=10)
    labeled = skimage.measure.label(cleared)
    props = skimage.measure.regionprops(labeled)
    pruned = []
    for j in im_thresh:
        subsample = j * cleared
        large = skimage.morphology.remove_small_objects(subsample, min_size=200)
        filled = scipy.ndimage.binary_fill_holes(large)
        pruned.append(filled)

    # Get the bounding boxes for objects < 100 µm^2. 
    bbox = [prop.bbox for prop in props]

    # Expand the bounding box by 2 px in each dimension
    crp = [np.s_[b[0]-3:b[2]+3, b[1]-3:b[3]+3] for b in bbox]

    # Set up the dataframe for final area storage.
    # Iterate through each microcolony and compute the area.
    for k, col in enumerate(crp):
        time = 0
        for j, t in enumerate(pruned):
            # Isolate the region. 
            _cleared = skimage.segmentation.clear_border(t[col])
            time += INTERVAL
            area = np.sum(_cleared)
            df = df.append({'date':DATE, 'carbon':CARBON, 'temp_C':TEMP,
                           'operator':OPERATOR, 'colony_idx':col_no, 
                           'area':area, 'time_min':time}, 
                           ignore_index=True) 
        col_no += 1

df.to_csv('output/{}_r{}_{}C_{}_{}_growth.csv'.format(DATE, RUN_NO, TEMP, CARBON,
                                                    OPERATOR), index=False)
conts = skimage.measure.find_contours(pruned[-1], level=0)
mwc.viz.growth_animation(im_filt, 
            'output/{}_r{}_{}C_{}_{}_growth_xy{}.gif'.format(DATE, RUN_NO, TEMP, CARBON, OPERATOR, i),
            contours=conts, descriptors={'bar_length':20, 'ip_dist':IP_DIST,
            'temperature':TEMP, 'carbon':CARBON, 'time_interval':INTERVAL})
