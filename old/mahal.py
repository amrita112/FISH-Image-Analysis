# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 14:09:18 2020
Mahalanobis distance from a point cloud in 3D


Input:
    img_dict: Dictionary in which each key is the name of a channel and 
              the corresponding value is the image of shape [n_planes, height, width] for that channel
    pixels_dict: Dictionary in which each key is the name of a group of pixels and 
                 the corresponding value is an array of pixel coordinates of shape [n_pixels, 3]
    origin_group: Name of the group of pixels to be used as the origin to calculate distance from (string)
    colors_dict: (optional) Dictionary in which each key is the name of a group of pixels and 
                 the corresponding value is a color
    
@author: Amrita Singh
"""


from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
from os.path import sep
from scipy.spatial.distance import mahalanobis

# %%
def mahal(img_dict, pixels_dict, origin_group, colors_dict = None, 
                     title = None, save = False, save_loc = None, save_file = None, log_scale = False):
    
    
    # %%
    channel_names = list(img_dict.keys())
    n_channels = len(channel_names)
    channels = range(n_channels)
    
    px_groups = list(pixels_dict.keys())
    
    
    # %% Get values for each channel for each pixel group
    px_vals = {}
    
    for group in px_groups:
        
        px = pixels_dict[group]
        
            
        px_vals[group] = {}
        
        for channel in channel_names:
            px_vals[group][channel] = img_dict[channel][(px[:, 0], px[:, 1], px[:, 2])]
    
            
    # %% Get mean and covariance matrix of origin group
    mu = np.zeros(n_channels)
    origin_px = np.zeros([n_channels, px_vals[origin_group][channel_names[0]].shape[0]])
    for c in range(n_channels):
        mu[c] = np.mean(px_vals[origin_group][channel_names[c]])
        origin_px[c, :] = px_vals[origin_group][channel_names[c]]

    cov = np.cov(origin_px) 
    VI = np.linalg.inv(cov)     
    
    
    # %%
    mahal = {}
    for group in px_groups:
        
        
        n_px = px_vals[group][channel_names[0]].shape[0]
        mahal[group] = np.zeros(n_px)
        u = np.zeros(n_channels)
        
        for px in range(n_px):
            
            for c in range(n_channels):
                u[c] = px_vals[group][channel_names[c]][px]
                
            mahal[group][px] = mahalanobis(u, mu, VI)
        
# %%
            
    if save:
        if save_loc == None:
            print('Please provide save location')
        else:
            if save_file == None:
                print('Please provide save file name')
            else:
                plt.savefig('{0}{1}{2}'.format(save_loc, sep, save_file))

















