# -*- coding: utf-8 -*-
# %%
import os
import sys
import numpy as np
import skimage.io
import glob
sys.path.insert(0, '../')
import mwc.image
import tqdm
root = '/Volumes/GDC_DATA_2/MBL_2018/20180618/'

# Set up the dictionary for the images.
iterators = [range(1, 11), range(11, 21), range(21, 31), range(31, 40)]
carbon_sources = ['glucose', 'RDM', 'acetate', 'glycerol']


for i in range(2,41):
    out = '{}/xy{:02d}/renamed'.format(root, i)
    if os.path.exists(out) == False:
        os.mkdir(out)
    
    files = glob.glob('{}xy{:02d}/*.tif'.format(root, i))
    files = np.sort(files)
    ims = skimage.io.ImageCollection(files)
    aligned = mwc.image.correct_drift(ims, crop=True)
    if i <=10:
        carbon = 'glucose'
    elif (i>10) & (i<=20):
        carbon = 'RDM'
    elif (i>20) & (i<=30):
        carbon = 'acetate'
    elif (i>30):
        carbon = 'glycerol'
    # Set the new base name.
    base_name = '{}/xy{:02d}/renamed/20180618_{}_xy{:02d}'.format(root, i, carbon, i)

    # Iterate through each image and rename.
    for k, im in enumerate(aligned):
        skimage.io.imsave('{}_t{:03d}.tif'.format(base_name, k), im)
       



