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
    bin_group: Name of the group of pixels to be binarized based on distance
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
def m_dist(img_dict, pixels_dict, origin_group, bin_group, thresh_scale = 2, colors_dict = None, 
                     title = None, save = False, save_loc = None, save_file = None, log_scale = False):
    
    
    # %%
    channel_names = list(img_dict.keys())
    n_channels = len(channel_names)
    
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
    m_dist = {}
    for group in px_groups:
        
        
        n_px = px_vals[group][channel_names[0]].shape[0]
        m_dist[group] = np.zeros(n_px)
        u = np.zeros(n_channels)
        
        for px in range(n_px):
            
            for c in range(n_channels):
                u[c] = px_vals[group][channel_names[c]][px]
                
            m_dist[group][px] = mahalanobis(u, mu, VI)
    
    # %% Set threshold
    
    all_points = m_dist[bin_group]
    center = np.median(all_points)
    neg_points = all_points[all_points < center]
    pos_points = neg_points + center
    all_points = np.concatenate([neg_points, pos_points])
    thresh = center + np.std(all_points)*thresh_scale
    
    # %%
    
    plt.figure()
    for group in px_groups:
        
        plt.hist(m_dist[group], 100, color = colors_dict[group], alpha = 0.5, label = group)
    
    y = np.linspace(plt.gca().get_ylim()[0], plt.gca().get_ylim()[1], 20)
    x = np.ones(20)*thresh
    plt.plot(x, y, linestyle = '--', color = 'r', label = 'Binarization threshold')
    
    plt.xlabel('Mahalanobis distance from {0} pixels'.format(origin_group), fontsize = 18)
    plt.ylabel('# pixels', fontsize = 18)
    plt.legend(fontsize = 18)
    plt.yscale('log')
    
# %%
            
    if save:
        if save_loc == None:
            print('Please provide save location')
        else:
            if save_file == None:
                print('Please provide save file name')
            else:
                plt.savefig('{0}{1}{2}'.format(save_loc, sep, save_file))


    return {'m_dist': m_dist, 'thresh': thresh, 'cov': cov}














