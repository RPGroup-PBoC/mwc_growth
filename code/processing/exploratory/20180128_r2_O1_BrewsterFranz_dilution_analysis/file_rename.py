import glob
import os
import shutil
import skimage.io
import pboc.image
import skimage.io
import numpy as np

# Define the experiment parameters.
DATE = 20121115
BASENAME = '37C_glucose_O1'
channels = ['p', 'y', 'r']
channel_dict = {c: i + 1 for i, c in enumerate(channels)}
keep_origs = True
generate_flatfield = True

# ----------------------------------------------------------------------------
# Nothing below here should change.
# ----------------------------------------------------------------------------
data_dir = '../../../data/images/{0}_{1}_dilution/'.format(DATE, BASENAME)
# Assemble the mean illumination image if appropriate.
flatfield = {}
field_ims = {}
if generate_flatfield is True:
    for ch in channel_dict.keys():
        if ch != 'p':
            # Load the slide images.
            slide_ims = glob.glob(
                '{0}{1}_fluorescent_slide*/*.tif'.format(data_dir, DATE))
            ims = skimage.io.ImageCollection(slide_ims)
            mean_im = pboc.image.projection(ims, mode='mean', median_filt=True)
            field_ims[ch] = mean_im
            flatfield[ch] = True
        else:
            flatfield[ch] = False
else:
    for ch in channel_dict.keys():
        flatfield[ch] = False


dark_avg = {}
for ch in channel_dict.keys():
    if ch != 'p':
        dark_files = glob.glob(
            '{0}{1}_camera_noise*/*.tif'.format(data_dir, DATE))
        dark_ims = skimage.io.ImageCollection(dark_files)
        dark_avg[ch] = pboc.image.projection(
            dark_ims, mode='mean', median_filt=True)

# %%
# Scrape positions.
samples = ['growth*/']
for s in samples:
    sub_samples = glob.glob('{0}{1}_{2}'.format(data_dir, DATE, s))
    for j, sub in enumerate(sub_samples):
        images = glob.glob('{0}/*.tif'.format(sub, s))
        for p in images:
            # Parse the position.
            name = p.split('/')[-1]
            if 'nofluo' in name:
                strain = 'autofluorescence'
                t = 0
                xy = int(name.split('_')[1].split('-')[0])
                ch = name.split('-')[1]

            else:
                strain = 'dilution'
                xy = int(name.split('_')[0].split('Pos')[1])
                ch = name.split('-')[1]
                t = int(name.split('-')[-1].split('.')[0])

            # Define the new name and copy.
            new_name = '{0}_{1}_t{2:05d}xy{3:03d}c{4}.tif'.format(
                DATE, BASENAME, t, xy, channel_dict[ch])
            if flatfield[ch] is True:
                im = skimage.io.imread(p)
                zero_im = np.zeros_like(im)
                ff_im = pboc.image.generate_flatfield(
                    im, dark_avg[ch], field_ims[ch], median_filt=True)
                skimage.io.imsave('{0}{1}'.format(sub, new_name), ff_im)
            else:
                shutil.copy(p, '{0}{1}'.format(sub, new_name))

    # # Move all original files into a separate folder.
    # if os.path.isdir('{0}/originals'.format(data_dir)) is False:
    #         os.mkdir('{0}/originals'.format(data_dir))
    #
    # shutil.move(p, '{0}'.format(sub, 'originals'))
    #
name
