import numpy as np
import skimage.io
import glob
import os
import shutil
skimage.io.use_plugin('freeimage')
data_dir = '../../../data/images/'
DATE = 20180123
BASE_NAME = 'artemis_37C_glucose_O2'
FOLDERS = ['{0}_snaps_deltaLacI_0ngmL_1'.format(
    DATE), '{0}_snaps_autofluorescence_0ngmL_1'.format(DATE),
    '{0}_growth_1'.format(DATE), '{0}_growth_fluo_1'.format(DATE)]
NUM_CHANNELS = [3, 3, 1, 2]
num_channels = {i: j for i, j in zip(FOLDERS, NUM_CHANNELS)}
snap_channel_names = {1: 'Brightfield', 2: 'YFP', 3: 'TRITC'}
growth_channel_names = {1: 'Brightfield', 2: 'TRITC'}
files = {i: j for i, j in zip(FOLDERS, NUM_CHANNELS)}

for i, f in enumerate(FOLDERS):
    tiffs = glob.glob('{0}{1}_{2}_dilution/{3}/*.tif'.format(data_dir, DATE,
                                                             BASE_NAME, f))

    for j, ome in enumerate(tiffs):
        # Find the positions.
        pos = int(ome.split('Pos')[1].split('.')[0])
        root_dir = '{0}{1}_{2}_dilution/{3}'.format(
            data_dir, DATE, BASE_NAME, f)
        if os.path.isdir('{0}/Pos{1}'.format(root_dir, pos)) == False:
            os.mkdir('{0}/Pos{1}'.format(root_dir, pos))
        if os.path.isdir('{0}/Pos{1}/ome_files'.format(root_dir, pos)) == False:
            os.mkdir('{0}/Pos{1}/ome_files'.format(root_dir, pos))
        # Load the image.
        im = skimage.io.imread(ome)
        im_shape = np.shape(im)

        # Determine the number of channels.
        n_chan = num_channels[f]

        if 'snaps' in f:
            time = 0
            for k in range(n_chan):
                # Save each image with the appropriate name.
                new_name = 'img_000000000_{0}_000.tif'.format(
                    snap_channel_names[k + 1])
                skimage.io.imsave('{0}/Pos{1}/{2}'.format(root_dir, pos,
                                                          new_name), im[:, :, k])
        elif 'growth_1' in f:
            for k in range(np.min(np.shape(im))):
                new_name = 'img_{0:09d}_Brightfield_000.tif'.format(k)
                skimage.io.imsave(
                    '{0}/Pos{1}/{2}'.format(root_dir, pos, new_name), im[k, :, :])

        elif 'growth_fluo' in f:
            for k in range(2):
                new_name = 'img_000000024_{0}_000.tif'.format(
                    growth_channel_names[k + 1])
                skimage.io.imsave(
                    '{0}/Pos{1}/{2}'.format(root_dir, pos, new_name), im[k, :, :])
        fname = ome.split('/')[-1]
        shutil.move(
            ome, '{0}/Pos{1}/ome_files/{2}'.format(root_dir, pos, fname))
