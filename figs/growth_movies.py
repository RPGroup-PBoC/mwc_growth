# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import skimage.io
sys.path.insert(0, '../')
import mwc.viz
colors = mwc.viz.personal_style()


fig, ax = plt.subplots(3, 3, sharex=True, sharey=True)
axs = ax.ravel()
for a in axs:
    a.set_xlim([-10, 500])
    a.set_ylim([-10, 500])
    a.set_xticks([])
    a.set_yticks([])
    a.set_facecolor('#000000')


# Define the ranges
carbon_sources = {'glycerol':(1, 2, 3), 'glucose':(11, 12, 13), 'RDM':(21, 23, 25), 'acetate':(34, )}

carbon_sources.items()

for carbon, pos in carbon_sources.items()
    

for i in range(3):
    files = glob.glob('/Volumes/GDC_DATA_2/MBL_2018/20180618/xy0{}/renamed/*.tif'.format(i+1))
    files = np.sort(files)
    ims = skimage.io.ImageCollection(files)
    ax[0, i].imshow(ims[-1], cmap=plt.cm.magma)

for i in range(3):
    files = glob.glob('/Volumes/GDC_DATA_2/MBL_2018/20180618/xy{:02d}/renamed/*.tif'.format(i+10))
    files = np.sort(files)
    ims = skimage.io.ImageCollection(files)
    ax[1, i].imshow(ims[-1], cmap=plt.cm.magma)

for i in range(3):
    files = glob.glob('/Volumes/GDC_DATA_2/MBL_2018/20180618/xy{:02d}/renamed/*.tif'.format(3+20))
    files = np.sort(files)
    ims = skimage.io.ImageCollection(files)
    ax[2, i].imshow(ims[-1], cmap=plt.cm.magma)

plt.subplots_adjust(hspace=-1, wspace=-1)
plt.tight_layout()
