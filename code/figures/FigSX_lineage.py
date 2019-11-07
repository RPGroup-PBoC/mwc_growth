#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import skimage.io
import mwc.viz
import mwc.io
import glob
import scipy.io
import mwc.stats
colors, _ = mwc.viz.personal_style()
#%%
# Load the phase contrast images and the segs
fluo_files = np.sort(glob.glob('../../data/20181021_r1_37C_glucose_O2_growth/xy4/fluor1/*.tif'))
seg_files = np.sort(glob.glob('../../data/20181021_r1_37C_glucose_O2_growth/xy4/seg/*_err.mat'))

# Load a representative seg file
seg1 = scipy.io.loadmat(seg_files[-1])
fluo1 = skimage.io.imread(fluo_files[-1]) 


s = np.s_[1250:1500, 0:400]
plt.imshow(seg1['mask_cell'] * fluo1)



# %%
