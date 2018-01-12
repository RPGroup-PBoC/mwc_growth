import numpy as np
import matplotlib.pyplot as plt
import tqdm
import pandas as pd
import skimage.io
import skimage.morphology
import skimage.measure
import skimage.filters
import scipy.ndimage
import statsmodels.tools.numdiff as smnd
import scipy.optimize


def compute_mean_bg(phase_image, fluo_image, method='isodata', obj_dark=True):
    """
    Computes the mean background fluorescence of the inverted segmentation
    mask.

    Parameters
    ----------
    phase_image : 2d-array, int or float.
        The phase contrast image used for generating the inverse segmentation
        mask. If this image is not a float with pixel values in (0, 1), it
        will be renormalized.
    fluo_image : 2d-array, int
        The fluorescence image used to calculate the mean pixel value. If
        flatfield correction is necessary, it should be done before this
        sending to this function.
    method: string, ['otsu', 'yen', 'li', 'isodata'], default 'isodata'
        Automated thresholding method to use. Default is 'isodata' method.
    obj_dark : bool, default True
        If True, objects will be **darker** than the automatically generated
        threshold value. If False, objects are deemed to be brighter.

    Returns
    -------
    mean_bg: float
        The mean background fluorescence of the image.
    """

    # Ensure that the image is renormalized.
    if (phase_image > 1.0).any():
        phase_image = (phase_image - phase_image.min()) /\
                      (phase_image.max() - phase_image.min())
    # Perform the background subtraction.
    im_blur = skimage.filters.gaussian(phase_image, sigma=50)
    im_sub = phase_image - im_blur

    # Determine the method to use.
    methods = {'otsu': skimage.filters.threshold_otsu,
               'yen': skimage.filters.threshold_yen,
               'li': skimage.filters.threshold_li,
               'isodata': skimage.filters.threshold_isodata}

    # Determine the threshold value.
    thresh_val = methods[method](im_sub)

    # Generate the inverted segmentation mask and dilate.
    if obj_dark is True:
        im_thresh = im_sub < thresh_val
    else:
        im_thresh = im_sub > thresh_val

    selem = skimage.morphology.disk(20)
    im_dil = skimage.morphology.dilation(im_thresh, selem=selem)

    # Mask onto the fluroescence image and compute the mean background value.
    mean_bg = np.mean(fluo_image[im_dil < 1])
    return mean_bg


def normalize_image(im, sub_bg=True):
    """
    Rescales the values of an image between 0 and 1. Can also perform a
    background subtraction.

    Parameters
    ----------
    im : 2d-array
        Image to be normalized.
    sub_bg: bool, default True.
        If True, a gaussian background subtraction is performed with
        a small sd.
    Returns
    -------
    im_norm : 2d-array
        Normalized image. If sub_bg is True, these values are on
        the domain [-1, 1]. If sub_bg is False, values are on [0, 1]
    """
    im_norm = (im - im.min()) / (im.max() - im.min())
    if sub_bg is True:
        im_blur = skimage.filters.gaussian(im_norm, sigma=5)
        im_norm = im_norm - im_blur
    return im_norm


def threshold_phase(im, min_int=0.15):
    """
    Performs an intensity based segmentation of a phase contrast image.
    This function uses Otsu's method to determine the threshold value.

    Parameters
    ----------
    im: 2d-array
        Image to be segmented. Desired objects in this image are assumed
        to be dark.
    min_int : float
        The maximum mean pixel intensity of a segmented object. This
        value must be between 0 and 1. Default is 0.15

    Returns
    -------
    mask: 2d-array, int
        Segmented image with labeled regions.
    """

    # Preprocess the phase image.
    im_sub = normalize_image(im)
    im_float = normalize_image(im, sub_bg=False)

    # Use Otsu's method.
    thresh = skimage.filters.threshold_otsu(im_sub)

    # Clean initial segmentation.
    seg = skimage.segmentation.clear_border(im_sub < thresh)
    seg = skimage.morphology.remove_small_objects(seg)
    mask = skimage.measure.label(seg)

    # Oversegment to correct for slight drift.
    selem = skimage.morphology.disk(2)
    mask = skimage.morphology.dilation(mask, selem)
    lab = skimage.measure.label(mask)

    # Impose minimum intensity filter.
    props = skimage.measure.regionprops(lab, im_float)
    final_im = np.zeros_like(mask)
    for prop in props:
        mean_int = prop.min_intensity
        if mean_int <= min_int:
            final_im += (lab == prop.label)
    mask = skimage.measure.label(final_im)
    return mask


def deterministic_log_post(alpha, I1, I2, p=0.5, neg=True):
    """
    Computes the log posterior for the deterministic measurement
    of the calibration factor α.

    Parmeters
    ---------
    alpha : float or int
        The calibration factor in units of A.U per molecule. This parameter
        must be positive.
    I1, I2 : arrays of floats or ints.
        The intensities of the two daugher cells. Each should be provided
        individually as 1d-arrays.
    p : float
        The probability of partitioning proteins into the first daughter cell.
        Defalt value is 0.5.
    neg : bool
        If True, the negative log posterior will be returned.

    Returns
    -------
    log_post: 1d-array
        The value of the log posterior at a given value of α.
    """
    # Ensure positivity of samples.
    if alpha < 0:
        return -np.inf
    if p < 0 or p > 1:
        raise RuntimeError('p must be between 0 and 1.')

    # Set the prior
    k = len(I1)
    prior = -k * np.log(alpha)

    # Determine the protein copy numbers deterministically.
    n1 = I1 / alpha
    n2 = I2 / alpha
    ntot = n1 + n2

    # Compute the binomial coefficient using log gamma functions.
    binom = scipy.special.gammaln(ntot + 1).sum() -\
        scipy.special.gammaln(n1 + 1).sum() - \
        scipy.special.gammaln(n2 + 1).sum()

    # Compute the probability factor of the binomial.
    prob = n1.sum() * np.log(p) + n2.sum() * np.log(1 - p)

    # Determine if the negative log posterior is desired.
    if neg is True:
        prefactor = -1
    else:
        prefactor = 1

    # Return the desired quantity.
    return prefactor * (prior + prob + binom)


def log_posterior_biexp(params, time, data, neg=True):
    """
    Computes the log posterior for a single exponential decay.

    Parameters
    ----------
    params : list or tuple
        The parameter values for the model. These should be in the
        order of beta, I_0_1, tau_1, I_0_2, tau_2.
    time: list or 1d-array
        The time points overwhich the bleaching occurs.
    data : list or 1d-array
        The bleaching data
    neg : bool
        If True, the negative log posterior is returned. Default is True.

    Returns
    -------
    log_post : float
        The value of the log posterior given the parameter inputs.
    """
    if neg is True:
        prefactor = -1
    else:
        prefactor = 1

    if (params < 0).any():
        return prefactor * -np.inf

    beta, I_0_1, tau_1, I_0_2, tau_2 = params

    k = len(data)
    mu = beta + I_0_1 * np.exp(-time / tau_1) + I_0_2 * np.exp(-time / tau_2)
    log_post = -np.log(tau_1) - np.log(tau_2) - (k / 2) * \
        np.log(np.sum((mu - data)**2))
    return prefactor * log_post


def estimate_calibration_factor(I1, I2, p=0.5):
    """
    Estimates the optimal value of α for a given data set by minimization.

    Parameters
    ----------
    I1, I2 : 1d-arrays
        The intensities of the two daughter cells. These should be provided
        individually as 1d-arrays.
    p : float
        The probability of paritioning into the first daughter cell. Default
        value is even, 0.5.

    Returns
    -------
    alpha_opt : float
        Best-fit value for the value of α in units of A.U. per molecule.
    sd : float
        The standard deviation of the the gaussian approximation of the
        posterior.
    """

    # Minimize the negative log posterior.
    popt = scipy.optimize.minimize_scalar(deterministic_log_post,
                                          args=(I1, I2, p, True))
    alpha_opt = popt.x

    # Compute the hessian.
    hess = smnd.approx_hess([alpha_opt], deterministic_log_post,
                            args=(I1, I2, p, False))
    cov = -np.linalg.inv(hess)
    sd = np.sqrt(cov[0])

    return [alpha_opt, sd[0]]


def cell_to_dict(file, eng, add_props=None, excluded_props=None):
    """
    Reads a single cell file and produces a dictionary containing
    the properties of interest.

    The returned properties are
    * birth - frame number at which the cell was born.
    * death - frame number at which the cell died.
    * divide - bool for an observed cell division.
    * ID - integer ID number of the cell.
    * motherID - integer ID number of the mother cell.
    * sisterID - integer ID number of the sister cell.
    * birth_fluo - fluorescence value at the cell's birth.
    * death_fluo - fluorescence value at the cell's death.
    * daughter_1_ID - integer ID number of the first daughter.
    * daughter_2_ID - integer ID number of the second daughter.


    Parameters
    ----------
    file: str
        Path of the cell file. This must be in a `.mat` format.
    eng: MATLAB engine object
        Engine of running matlab session.
    add_props : dict, default None
        Dictionary of additional properties (not found in the mat file)
        to be included in the returned dictionary.
    excluded_props: list of str
        Properties of cell.mat file to be ignored. These must be
        exactly how they are defined in the cell file.

    Returns
    -------
    cell_dict : dictionary
        Dictionary of all extracted properties from the cell files.
    """

    # Ensure the supplied file is actually a .mat and other types are correct.
    if file.split('.')[-1] != 'mat':
        raise TypeError("supplied file {0} is not a `.mat` file.".format(file))
    if add_props is not None and type(add_props) is not dict:
        raise TypeError(
            "add_props is {0} and not dict.".format(type(add_props)))
    if excluded_props is not None and type(excluded_props) is not list:
        raise TypeError(
            "add_props must be list. Type is currently {0}.".format(type(excluded_props)))

    # Define the values of interest.
    vals = ['birth', 'death', 'divide', 'ID', 'motherID', 'sisterID',
            'daughter_1_ID', 'daughter_2_ID', 'birth_fluo', 'death_fluo', 'birth_area', 'death_area']

    # Load the mat file using MATLAB.
    eng.workspace['f'] = file
    mat = eng.eval('load(f)')

    # Assemble the dictionary for constant properties.
    cell_dict = {v: mat[v] for v in vals[:-6]}
    daughters = np.array(mat['daughterID'])

    # Determine  if daughters were produced. If not, change ID to NaN.
    if len(daughters) == 0:
        daughter_1, daughter_2 = None,  None
    else:
        daughter_1, daughter_2 = daughters[0]
    cell_dict['daughter_1_ID'] = daughter_1
    cell_dict['daughter_2_ID'] = daughter_2

    # Extract fluorescence information -- This is a bit gross but checked.
    # Get number of fluorescence channels.
    fluo_channels = [f for f in mat['CellA'][0].keys() if 'fl' in f]
    n_channels = int(len(fluo_channels) / 3)
    for n in range(n_channels):
        _n = n + 1
        try:
            fluo = [mat['CellA'][i]['fl{0}'.format(
                _n)]['sum'] for i, _ in enumerate(mat['CellA'])]
            nonzero = np.where(np.array(fluo) > 0)[0]
            num_exposures = len(nonzero)
            cell_dict['fluor{0}_birth_fluo'.format(_n)] = fluo[nonzero.min()]
            cell_dict['fluor{0}_death_fluo'.format(_n)] = fluo[nonzero.max()]
            cell_dict['birth_area'] = mat['CellA'][nonzero.min()]['coord']['A']
            cell_dict['death_area'] = mat['CellA'][nonzero.max()]['coord']['A']
        except:
            cell_dict['fluor{0}_birth_fluo'.format(_n)] = 0
            cell_dict['fluor{0}_death_fluo'.format(_n)] = 0
            cell_dict['birth_area'] = mat['CellA'][0]['coord']['A']
            cell_dict['death_area'] = mat['CellA'][-1]['coord']['A']
            num_exposures = 0
            cell_dict['fluor{0}_num_exposures'.format(_n)] = num_exposures

    # Deal with exclusion and addition of props.
    if excluded_props is not None:
        new_dict = {}
        keys = cell_dict.keys()
        for key in keys:
            if key not in excluded_props:
                new_dict[key] = cell_dict[key]
        cell_dict = new_dict
    if add_props is not None:
        for key in add_props.keys():
            cell_dict[key] = add_props[key]

    # Return the cell dictionary.
    return cell_dict


def parse_cell_files(files, eng, verbose=False, **kwargs):
    """
    Executes cell_to_dict across a list of files and returns a Pandas DataFrame.
    """
    if type(files) is not list:
        raise TypeError("'files' is type {0} not list.".format(type(files)))
    if verbose:
        files = tqdm.tqdm(files)
    for i, f in enumerate(files):
        cell_dict = cell_to_dict(f, eng, **kwargs)
        if i == 0:
            keys = cell_dict.keys()
            df = pd.DataFrame([], columns=keys)
            df = df.append(cell_dict, ignore_index=True)
        else:
            df = df.append(cell_dict, ignore_index=True)
    return df
