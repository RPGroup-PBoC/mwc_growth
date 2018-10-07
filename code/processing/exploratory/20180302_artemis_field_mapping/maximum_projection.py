import numpy as np
import skimage.io
import sys
import tqdm
import glob

# Load the images. 
files = glob.glob('../../../../data/images/20180302_artemis_spatial_profiling/single_bead_1/*/*.tif')
print(len(files))
for i, f in enumerate(tqdm.tqdm(files)):
    if i == 0:
        proj = skimage.io.imread(f)
    else:
       _im = skimage.io.imread(f)
       proj = np.max(np.array([_im, proj]), axis=0)
skimage.io.imsave('output/max_projection.tif', proj)
