# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 15:11:16 2020

Filter a 2D image with a 2D difference of gaussians filter

@author: Amrita Singh
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from skimage.filters import gaussian

from scipy.stats import norm

def diff_gauss(sigma_small, sigma_large, img, plot_kernel = 0):

    #sigma_large = 10 # In pixels
    #sigma_small = 3 # In pixels

    #x = np.linspace(norm.ppf(0.01), norm.ppf(0.99), 1000)
    x = np.linspace(-sigma_large*5, sigma_large*5, 1000)
    small = norm.pdf(x, scale = sigma_small)
    large = norm.pdf(x, scale = sigma_large)
    diff = small - large

    im_diff = np.zeros(img.shape)

    im_filt_large = gaussian(img, sigma = sigma_large)
    im_filt_small = gaussian(img, sigma = sigma_small)
    im_diff = im_filt_small - im_filt_large

    if plot_kernel:

        plt.figure()
        plt.plot(x, diff, label = 'diff')
        plt.plot(x, small, label = 'Small sigma filter')
        plt.plot(x, large, label = 'Large sigma filter')
        plt.xlabel('Distance from center of kernel (px)')
        plt.ylabel('Kernel')
        plt.title('Difference of gaussians filter')
        plt.legend()

    return im_diff
