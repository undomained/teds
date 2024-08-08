#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 2023.

@author: Manu Goudar


"""

import numpy as np
from skimage.measure import label


class PlumeCont:
    def __init__(self, flag, segmentedimage, plume):
        self.flag = flag
        self.segmentimage = segmentedimage
        self.plume = plume


def cdfthresholdbasedplume(img, cdfweight):
    """_summary_

    Args:
        img (matrix, float): Image
        cdfweight (float): Weight in Cumulative Distributve Function.
    Returns:
        plume: Plume image
    """
    # compute a histogram to define threshold
    hist = np.histogram(img.ravel(), bins=100)
    cdf = np.cumsum(hist[0])
    cdf = cdf / cdf[-1]
    # Value of threshold
    threshold = np.interp(cdfweight, cdf, (hist[1][1:] + hist[1][:-1]) * 0.5)
    # plume image
    return img >= threshold


def _extractuniquelabelsaroundsource(orig_id, segmented_img, nopixels=2):
    # Extract all labels around the source
    i1 = orig_id[0] - nopixels
    i2 = orig_id[0] + nopixels + 1
    j1 = orig_id[1] - nopixels
    j2 = orig_id[1] + nopixels + 1
    return np.unique(segmented_img[i1:i2, j1:j2])


def _extractlargestplume(label_ids, segmented_img):
    maxpixels = 0
    maxpixels_id = 0
    for i in label_ids:
        if i != 0:  # if it is not background pixels
            tmp = (segmented_img == i).sum()  # get number of pixels
            if tmp > maxpixels:
                maxpixels = tmp
                maxpixels_id = i
    return maxpixels, maxpixels_id


def get_thresholdbasedplume(img, orig_id, cdfweight, minimumplumepixels=7):
    # compute plume based on a threshold
    plume_img = cdfthresholdbasedplume(img, cdfweight)

    # Label computed plume
    segmented_img = label(plume_img, background=False)

    # Extract the plume
    label_ids = _extractuniquelabelsaroundsource(orig_id, segmented_img, nopixels=2)

    # here based on number of labels around origin define the plume.
    if (len(label_ids) == 1) & (label_ids[0] == 0): 
        # There is only one label which is background.
        return PlumeCont(False, segmented_img, plume_img)
    else:  # There are more labels around origin and one of them can be plume
        maxpixels, maxpixels_id = _extractlargestplume(label_ids, segmented_img)
        # If the number of pixels are more than the minimum number of plume pixels
        if maxpixels >= minimumplumepixels:
            # get its label
            finaldetectedplume = np.zeros_like(segmented_img, dtype=np.bool_)  # create a new label image
            finaldetectedplume[segmented_img == maxpixels_id] = True  # create a plume image
            return PlumeCont(True, segmented_img, finaldetectedplume)
        else:
            return PlumeCont(False, segmented_img, plume_img)
