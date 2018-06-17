import numpy as np
import skimage.io
import skimage.morphology
import skimage.filters
import skimage.segmentation
import skimage.feature
import scipy.ndimage
import pandas as pd
from mwc.image import projection, generate_flatfield, correct_drift
import pytest


# Set up sample arrays
ones_im = np.ones((5, 5)).astype('uint16')
threes_im = (np.ones((5, 5)) * 3).astype('uint16')
fives_im = (np.ones((5, 5)) * 5).astype('uint16')
tens_im = (np.ones((5, 5)) * 10).astype('uint16')
field_im = skimage.morphology.disk(2) + 1
flat_im = ((tens_im - ones_im) * np.mean(field_im - ones_im)) / \
    (field_im - ones_im)


def test_projection():
    im_array = [threes_im, fives_im, tens_im]
    assert (fives_im + 1 == projection(im_array, mode='mean',
                                       median_filt=False)).all()
    assert (fives_im == projection(im_array, mode='median',
                                   median_filt=False)).all()
    assert (threes_im == projection(im_array, mode='min',
                                    median_filt=False)).all()
    assert (tens_im == projection(im_array, mode='max',
                                  median_filt=False)).all()


def test_generate_flatfield():
    assert (flat_im == generate_flatfield(tens_im, ones_im, field_im,
                                          median_filt=False)).all()


def test_correct_drift():
    shift = (-1, 1)
    shifted_ones_im = scipy.ndimage.fourier_shift(ones_im, shift) 
    ones_im_list = [ones_im, shifted_ones_im]
    aligned_ones_ims = [np.append(ones_im, np.zeros((1, 5)), axis=0)[1:,:], shifted_ones_im]
    assert (correct_drift(ones_im_list) == aligned_ones_ims).all()
