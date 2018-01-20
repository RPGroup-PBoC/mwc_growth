import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '../../../')
import mwc_growth as mwc
import pboc.plotting
import pboc.image
import pandas as pd
import glob
import skimage.io
import scipy.optimize
from tqdm import tqdm
colors = pboc.plotting.set_plotting_style()
import imp
imp.reload(mwc)
# Define the experimental parameters.
DATE = 20180110
BASE_STRAIN = '37C_glucose_O2'
SAMPLE_NAMES = ['auto', 'delta']
root_dir = '../../../data/images/{0}_{1}_dilution/'.format(
    DATE, BASE_STRAIN)
IP_DIST = 0.16
EXPOSURE_MS = 50  # Exposure in ms.

# %% Load the flatfield correction.
slide_files = glob.glob(
    '{0}/*_TRITC_fluorescent_slide*/Pos*/*.tif'.format(root_dir))
slide_ims = skimage.io.ImageCollection(slide_files)
mean_field = np.mean(slide_ims, axis=0)

# %% Process the bleaching files.
col_names = ['strain', 'cell_id', 'area', 'total_intensity',
             'elapsed_time_s', 'mean_bg', 'date']
bleaching_df = pd.DataFrame([], columns=col_names)
print('Beginning processing...')
for i, nom in enumerate(tqdm(SAMPLE_NAMES)):
    files = glob.glob(
        '{0}photobleaching/*{1}*'.format(root_dir, nom))
    cell_counter = 0
    for j, pos in enumerate(tqdm(files, desc='Processing positions')):
        # Load the images
        try:
            file = glob.glob('{0}/*.tif'.format(pos))

            ims = skimage.io.ImageCollection(file)

            # Split in to phase and intensity.
            bf_im = ims[0][:, 0, :, :][0]
            fluo_ims = ims[0][:, 1, :, :]
        except:
            bf_file = glob.glob('{0}/Pos*/*Brightfield*.tif'.format(pos))
            fl_files = glob.glob('{0}/Pos*/*TRITC*.tif'.format(pos))
            bf_im = skimage.io.imread(bf_file[0])
            fluo_ims = skimage.io.ImageCollection(fl_files)

        # Segment the image and generate the inverse mask.
        mask = mwc.threshold_phase(bf_im)
        inv_mask = (mask < 1)

        # Flatten the fluorescence images.
        zero_im = np.zeros_like(bf_im)
        fluo_flat = [pboc.image.generate_flatfield(
            im, zero_im, mean_field) for im in fluo_ims]
        for k, im in enumerate(fluo_flat):
            print('Processing image {0} of {1}.'.format(
                k + 1, len(fluo_flat) + 1))
            # Compute the mean background of the position.
            mean_bg = np.mean(im[inv_mask])
            props = skimage.measure.regionprops(mask, im)
            cell_no = 0
            for z, p in enumerate(props):
                area = p.area * IP_DIST**2
                total_intensity = p.mean_intensity * area
                label = z + cell_counter
                elapsed_time = k * EXPOSURE_MS / 1E3
                cell_dict = dict(strain=nom, cell_id=label,
                                 total_intensity=total_intensity - area * mean_bg,
                                 elapsed_time_s=elapsed_time,
                                 mean_bg=mean_bg, date=DATE, area=area)
                bleaching_df = bleaching_df.append(
                    cell_dict, ignore_index=True)
        cell_counter += np.max(mask)
print('...bleaching files processed, begining analysis...')


# %% Rescale data
grouped = bleaching_df.groupby(['strain', 'cell_id'])
dfs = []
for g, d in grouped:
    i0 = d[d['elapsed_time_s'] == 0.0]['total_intensity'].values[0]
    d.loc[:, 'rescaled_intensity'] = d.loc[:, 'total_intensity'] / i0
    dfs.append(d)
rescaled_df = pd.concat(dfs, axis=0)
rescaled_df.to_csv('output/{0}_sfGFP_bleaching.csv'.format(DATE))

# %% Compute the mean traces for each.
grouped = rescaled_df.groupby(['elapsed_time_s'])
mean_df = pd.DataFrame([], columns=['mean_auto', 'mean_delta',
                                    'elapsed_time', 'sub_mean'])
for g, d in grouped:
    mean_auto = d[d['strain'] == 'autofluorescence']['total_intensity'].mean()
    mean_delta = d[d['strain'] == 'delta']['total_intensity'].mean()
    sub = mean_delta - mean_auto
    cell_dict = dict(elapsed_time=g, mean_auto=mean_auto,
                     mean_delta=mean_delta, sub_mean=sub)
    mean_df = mean_df.append(cell_dict, ignore_index=True)

mean_df['sub_rescaled'] = mean_df['sub_mean'] / mean_df.iloc[0]['sub_mean']

# %% Perform inference of autofluorescence bleaching to set informative priors.
guess = [0, 0.5, 10, 0.5, 20]
popt = scipy.optimize.minimize(mwc.log_posterior_biexp, guess, method='Powell',
                               args=(mean_df['elapsed_time'], mean_df['sub_rescaled'], True))

# Extract the parameters.
bg, beta_1, tau_1, beta_2, tau_2 = popt.x
fit_stats = pd.DataFrame(popt.x).T
fit_stats.columns = ['bg', 'beta_1', 'tau_1', 'beta_2', 'tau_2']

fit_stats.to_csv('output/{0}_{1}_bleaching_constants.csv'.format(DATE, BASE_STRAIN),
                 index=False)
# Plot the fit.
time_range = np.linspace(0, 100, 1000)
fit = bg + beta_1 * np.exp(-time_range / tau_1) + \
    beta_2 * np.exp(-time_range / tau_2)
fig, ax = plt.subplots(1, 1)
ax.set_xlabel('time [s]')
ax.set_ylabel('rescaled intensity')
_ = ax.plot(mean_df['elapsed_time'], mean_df['sub_rescaled'], '.', color='slategray',
            label='auto subtracted mean')
_ = ax.plot(time_range, fit, color='dodgerblue', label='biexponential fit')
plt.legend()
plt.tight_layout()
plt.savefig('output/{0}_bleaching_fit.png'.format(DATE))
