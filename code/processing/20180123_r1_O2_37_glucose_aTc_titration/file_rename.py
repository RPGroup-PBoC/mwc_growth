import glob
import os
import shutil
import skimage.io
import pboc.image
import skimage.io
import numpy as np

# Define the experiment parameters.
DATE = 20180123
BASENAME = '37C_glucose_O2'
channels = ['Brightfield', 'TRITC', 'YFP']
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
        if ch != 'Brightfield':
            # Load the slide images.
            slide_ims = glob.glob(
                '{0}{1}_{2}_fluorescent_slide*/Pos*/*.tif'.format(data_dir, DATE, ch))
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
    if ch != 'Brightfield':
        dark_files = glob.glob(
            '{0}{1}_{2}_camera_noise*/Pos*/*.tif'.format(data_dir, DATE, ch))
        dark_ims = skimage.io.ImageCollection(dark_files)
        dark_avg[ch] = pboc.image.projection(
            dark_ims, mode='mean', median_filt=True)

# %%
# Scrape positions.
samples = ['*snaps*/']
for s in samples:
    if 'snaps' in s:
        sub_samples = glob.glob('{0}{1}_{2}'.format(data_dir, DATE, s))
    else:
        sub_samples = glob.glob('{0}{1}_{2}'.format(data_dir, DATE, s))
    for j, sub in enumerate(sub_samples):
        positions = glob.glob('{0}/Pos*/'.format(sub, s))

        for p in positions:
            # Parse the position.
            folder = p.split('Pos')[0]
            xy = '{0:03d}'.format(int(p.split('/')[-2].split('s')[1]))

            # Loop through provided channels and rename images.
            for ch in channels:
                files = glob.glob(p + '*{0}*.tif'.format(ch))
                files.sort()
                for f in files:
                    # Parse the relevant info
                    _, t, c, _ = f.split('/')[-1].split('_')

                    # Define the new name and copy.
                    new_name = '{0}_{1}_t{2:05d}xy{3:03d}c{4}.tif'.format(
                        DATE, BASENAME, int(t), int(xy), channel_dict[ch])
                    if flatfield[ch] is True:
                        im = skimage.io.imread(f)
                        zero_im = np.zeros_like(im)
                        ff_im = pboc.image.generate_flatfield(
                            im, dark_avg[ch], field_ims[ch], median_filt=True)
                        skimage.io.imsave('{0}{1}'.format(
                            folder, new_name), ff_im)
                    else:
                        shutil.copy(f, '{0}{1}'.format(folder, new_name))

        # Rename and save the metadata files.
        mfiles = glob.glob('{0}metadata.txt'.format(p))
        if len(mfiles) != 0:
            for m in mfiles:
                new_m_name = 'xy{0:03d}_metadata.txt'.format(int(xy))
                shutil.copy(m, '{0}{1}'.format(data_dir, new_m_name))

        # Move all original files into a separate folder.
        if os.path.isdir('{0}/originals'.format(data_dir)) is False:
            os.mkdir('{0}/originals'.format(data_dir))

        for p in positions:
            shutil.move(p, '{0}/{1}/{2}'.format(sub,
                                                'originals', p.split('/')[-3]))

    # Determine if originals should be kept.
    if keep_origs is False:
        shutil.rmtree('{0}{1}'.format(data_dir, 'originals'))
