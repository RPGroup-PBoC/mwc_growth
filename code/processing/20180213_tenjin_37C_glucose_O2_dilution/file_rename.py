# -*- coding: utf-8 -*-
import numpy as np
import glob
import skimage.io
import skimage.morphology
import scipy.ndimage
import os
import shutil
import sys
import warnings
import tqdm
sys.path.insert(0, '../../../')
import mwc.image
skimage.io.use_plugin('freeimage')

# Set the experiment parameters.
DATE = 20180213
TEMP = 37  # in °C
CARBON = 'glucose'
MICROSCOPE = 'tenjin'
OPERATOR = 'O2'

# Define the channels of interest
channels = ['Brightfield', 'mCherry', 'YFP']
channel_dict = {c: i + 1 for i, c in enumerate(channels)}

# ################
# Nothing below should change from experiment to experiment.
# ################
data_dir = '../../../data/images/{0}_{1}_{2}C_{3}_{4}_dilution/'.format(
    DATE, MICROSCOPE, TEMP, CARBON, OPERATOR)
selem = skimage.morphology.square(3)

#%% Create flatfield images.
field_avgs = {}
noise_avgs = {}
ff_channels = ['mCherry', 'YFP']
ff_dict = {i + 2: ch for i, ch in enumerate(ff_channels)}
for ch in ff_channels:
    # Grab all of the ff images.
    noise_files = glob.glob('{0}{1}_camera_noise_w{2}*.TIF'.format(data_dir, DATE,
                                                                   channel_dict[ch] - 1))
    field_files = glob.glob('{0}{1}_fluorescent_slide_*{2}*.TIF'.format(data_dir, DATE,
                                                                        ch))
    # Load the images.
    noise_ims = skimage.io.ImageCollection(noise_files, conserve_memory=False)
    field_ims = skimage.io.ImageCollection(field_files, conserve_memory=False)

    # Create a mean projection of both.
    noise_avg = mwc.image.projection(noise_ims, mode='mean', median_filt=False)
    field_avg = mwc.image.projection(field_ims, mode='mean', median_filt=False)

    # Add them to the dictionary.
    noise_avgs[channel_dict[ch]] = noise_avg
    field_avgs[channel_dict[ch]] = field_avg


# %% Rename files to SuperSegger requirements.
samples = ['snaps', 'growth']
snap_groups = []
for i, s in enumerate(tqdm.tqdm(samples)):
    samp_files = glob.glob('{0}*{1}*.TIF'.format(data_dir, s))
    samp_files
    # Make sample folders if necessary.
    if os.path.isdir('{0}{1}'.format(data_dir, s)) == False:
        os.mkdir('{0}{1}'.format(data_dir, s))

    for j, f in enumerate(tqdm.tqdm(samp_files)):
        # Check if there are multiple wavelengths.
        multi_wave = np.sum(
            [i in f for i in channels]).astype('bool')

        # Get the particulars of the file.
        if multi_wave == True:
            if 'ngmL' in f:
                _, ident, strain, atc_conc, ch, pos = f.split(
                    '/')[-1].split('_')
                atc_conc = int(atc_conc.split('ngmL')[0])
            elif 'growth_fluo' in f:
                 _, ident, strain, ch, pos, _ = f.split('/')[-1].split('_')
                 atc_conc = 'mixed'
            else:
                _, ident, strain, ch, pos = f.split('/')[-1].split('_')
                atc_conc = 'mixed'
            ch = int(ch[1])
            pos = int(pos.split('s')[1].split('.TIF')[0])
            time = 0

        else:
            _, ident, pos, time = f.split('/')[-1].split('_')
            pos = int(pos.split('s')[1])
            time = int(time.split('t')[1].split('.TIF')[0])
            ch = 1
            strain = 'dilution'
            atc_conc = 'mixed'

        if s == 'snaps':
            snap_group = 'snaps_{0}_{1}ngmL'.format(strain, atc_conc)
            if snap_group not in snap_groups:
                snap_groups.append(snap_group)

        # Define the new name that's compatible with SuperSegger
        new_name = '{0}_{1}_{2}_{3}ngmL_t{4:05d}xy{5:03d}c{6}.tif'.format(
            DATE, ident, strain, atc_conc, time, pos, ch)

        # Peform the flatfield correction if necessary.
        if ch in ff_dict.keys():
            im = skimage.io.imread(f)
            ff_im = mwc.image.generate_flatfield(im, noise_avgs[ch], field_avgs[ch],
                                                 median_filt=False)
            ff_filt = scipy.ndimage.median_filter(ff_im, footprint=selem)

            # Convert to 16 bit.
            ff_im = np.round(ff_filt).astype(np.uint16)

            # Save image in correct folder.
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                skimage.io.imsave(
                    '{0}{1}/{2}'.format(data_dir, s, new_name), ff_im)
        else:
            shutil.copy(f, '{0}{1}/{2}'.format(data_dir, s, new_name))


# Make a directory for the originals and move them there.
if os.path.isdir('{0}originals'.format(data_dir)) == False:
    os.mkdir('{0}originals'.format(data_dir))
files = glob.glob('{0}*.TIF'.format(data_dir))
for f in files:
    shutil.move(f, '{0}originals'.format(data_dir))

# %% Adjust the growth_fluo timepoints to the correct time.
growth_files = glob.glob('{0}growth/*dilution*.tif'.format(data_dir))
growth_fluo_files = glob.glob('{0}growth/*growth_fluo*.tif'.format(data_dir))

# Find the maximum time point for the growth experiment.
max_time = int(np.sort(growth_files)
               [-1].split('/')[-1].split('_t')[-1].split('xy')[0])

# Iterate through each growth fluo file and change the time and strain.
for i, f in enumerate(growth_fluo_files):
    _, ident, _, atc, seqinfo = f.split('/')[-1].split('_')
    pos = seqinfo.split('xy')[1]
    new_name = '_'.join([str(DATE), ident, 'dilution', atc,
                         't{0:05d}xy{1}'.format(max_time, pos)])

    # Rename the file
    os.rename(f, '{0}growth/{1}'.format(data_dir, new_name))

#%% Complete the fluorescence sequence.
missing_channels = [2]
im = skimage.io.imread(growth_files[0])
zeros_im = np.zeros_like(im)
for i, f in enumerate(growth_files):
    # Get the name of the file.
    name = f.split('/')[-1]
    prefix = name[:-5]
    for ch in missing_channels:
        new_name = '{0}growth/{1}{2}.tif'.format(data_dir, prefix, ch)
        if os.path.exists(new_name) == False:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                skimage.io.imsave(new_name, zeros_im)

# %% Partition growth groups into superseggable chunks.
max_pos = int(np.sort(growth_files)
              [-1].split('/')[-1].split('xy')[1].split('c')[0])
num_groups = int(np.ceil(max_pos / 10))
if num_groups > 1:
    for i in range(1, num_groups):
        if os.path.isdir('{0}growth_{1}'.format(data_dir, i)) == False:
            os.mkdir('{0}growth_{1}'.format(data_dir, i))
        # Grab all of the files for the position chunks.
        files = glob.glob('{0}growth/*xy0{1}*'.format(data_dir, i))
        for f in files:
            shutil.move(f, '{0}growth_{1}/'.format(data_dir, i))

    # Rename the first growth folder.
    shutil.move('{0}growth'.format(data_dir), '{0}growth_0'.format(data_dir))

# %% Make and move snap group folders.
for i, s in enumerate(snap_groups):
    if os.path.isdir('{0}{1}'.format(data_dir, s)) == False:
        os.mkdir('{0}{1}'.format(data_dir, s))
    files = glob.glob('{0}snaps/*{1}*.tif'.format(data_dir, s))
    for f in files:
        shutil.move(f, '{0}{1}'.format(data_dir, s))

# Remove the snaps directory.
os.rmdir('{0}snaps'.format(data_dir))

print("""
Finished! Thank you come again.
      """)
