# -*- coding: utf-8 -*-
import skimage.segmentation
import skimage.io
import skimage.morphology
import skimage.filters
import skimage.measure
import skimage.feature
import scipy.ndimage
import numpy as np
import pandas as pd
import bokeh.io
import tqdm



def projection(im, mode='mean', median_filt=True):
    R"""
    Computes an average image from a provided array of images.

    Parameters
    ----------
    im : list or arrays of 2d-arrays
        Stack of images to be filtered.
    mode : string ('mean', 'median', 'min', 'max')
        Type of elementwise projection.
    median_filt : bool
        If True, each image will be median filtered before averaging.
        Median filtering is performed using a 3x3 square structural element.

    Returns
    -------
    im_avg : 2d-array
        Projected image with a type of int.
    """
    # Determine if the images should be median filtered.
    if median_filt is True:
        selem = skimage.morphology.square(3)
        im_filt = [scipy.ndimage.median_filter(i, footprint=selem) for i in im]
        im = im_filt
    # Get the image type
    im_type = im[0].dtype
    # Determine and perform the projection.
    if mode is 'mean':
        im_proj = np.mean(im, axis=0)
    elif mode is 'median':
        im_proj = np.median(im, axis=0)
    elif mode is 'min':
        im_proj = np.min(im, axis=0)
    elif mode is 'max':
        im_proj = np.max(im, axis=0)
    return im_proj.astype(im_type)


def generate_flatfield(im, im_dark, im_field, median_filt=True):
    """
    Corrects illumination of a given image using a dark image and an image of
    the flat illumination.

    Parameters
    ----------
    im : 2d-array
        Image to be flattened.
    im_field: 2d-array
        Average image of fluorescence illumination.
    median_filt : bool
        If True, the image to be corrected will be median filtered with a
        3x3 square structural element.

    Returns
    -------
    im_flat : 2d-array
        Image corrected for uneven fluorescence illumination. This is performed
        as

        im_flat = ((im - im_dark) / (im_field - im_dark)) *
                   mean(im_field - im_dark)

    Raises
    ------
    RuntimeError
        Thrown if bright image and dark image are approximately equal. This
        will result in a division by zero.
    """

    # Compute the mean field image.
    mean_diff = np.mean(im_field - im_dark)

    if median_filt is True:
        selem = skimage.morphology.square(3)
        im_filt = scipy.ndimage.median_filter(im, footprint=selem)

    else:
        im_filt = im

    # Compute and return the flattened image.
    im_flat = (im_filt - im_dark) / (im_field - im_dark) * mean_diff
    return im_flat

def correct_drift(images, verbose=True, crop=True):
    R"""
    Corrects xy drift between a set of images. 

    Parameters
    ----------
    images: list, array, or ImageCollection object
        List of images that are to be corrected. 
    verbose: bool
        If True, a progress par will be printed to screen during reistration.
    crop: bool
        If True, the interpolated regions will be cropped from the image.
    
    Returns
    -------
    shifted: list
        List of images with corrected positioning. Shifted images will be
        larger than those supplied
    """
    shifted = []
    if verbose == True:
        iterator = tqdm.tqdm(range(1, len(images)), desc='performing image regsitration...')
    else:
        iterator = range(1, len(images))

    # Set up counters for the maximjum magnitude shift.
    max_x = 0
    max_y = 0
    for i in iterator:

        # Complete the registration
        drift, _, _ = skimage.feature.register_translation(images[0], images[i])
        shift = scipy.ndimage.interpolation.shift(images[i], drift, mode='nearest')

        # Store the maximum magnitude shift
        y, x = drift
        if np.abs(x) > np.abs(max_x):
            max_x = x
        if np.abs(y) > np.abs(max_y):
            max_y = y
        shifted.append(shift)
    max_x = int(max_x)
    max_y = int(max_y)
    # Determine if the shifted images should be cropped
    if crop == True:
        if  (max_x < 0) & (max_y < 0):
            restriction = np.s_[0:max_y+1, 0:max_x+1]
        elif  (max_x < 0) & (max_y >= 0):
            restriction = np.s_[max_y:-1, 0:max_x+1]
        elif (max_x >=0 ) & (max_y < 0):
            restriction = np.s_[0:max_y+1, max_x:-1]
        elif (max_x >=0) & (max_y >=0):
            restriction = np.s_[max_y:-1, max_x:-1]

        _shifted = [im[restriction] for im in shifted]
        shifted = _shifted
    return shifted 

