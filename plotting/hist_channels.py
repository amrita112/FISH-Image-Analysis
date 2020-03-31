# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 08:07:47 2020

Input:
    img_dict: Dictionary in which each key is the name of a channel and 
              the corresponding value is the image of shape [n_planes, height, width] for that channel
    pixels_dict: Dictionary in which each key is the name of a group of pixels and 
                 the corresponding value is an array of pixel coordinates of shape [n_pixels, 3]
    colors_dict: (optional) Dictionary in which each key is the name of a group of pixels and 
                 the corresponding value is a color
    max_per_group: Default 10000. Maximum number of pixels per group to plot. 
                   If number of pixels exceeds this, pixels will be randomly sub-sampled.

@author: Amrita Singh
"""

from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
from os.path import sep

def hist_channels(img_dict, pixels_dict, max_per_group = 10000, bins = 100, colors_dict = None, 
                  do_bin = False, bin_group = None, thresh_scale = None, title = None, 
                  save = False, save_loc = None, save_file = None, log_scale = True):
    
    # %%
    channel_names = list(img_dict.keys())
    n_channels = len(channel_names)
    thresh = {}
    px_groups = list(pixels_dict.keys())
    
    
    # %% Get values for each channel for each pixel group
    px_vals = {}
    
    for group in px_groups:
        
        px = pixels_dict[group]
        n_px = px.shape[0]
        
        # Sub-sample if too many pixels to plot
        if n_px > max_per_group:
            idx = np.random.randint(int(n_px/max_per_group), size = n_px)
            px = px[idx == 0, :]
            
        px_vals[group] = {}
        
        for channel in channel_names:
            px_vals[group][channel] = img_dict[channel][(px[:, 0], px[:, 1], px[:, 2])]
            
    # %%
    fig, ax = plt.subplots(nrows = n_channels, ncols = 1, figsize = (10, 15))
    for c in range(n_channels):
        
        channel_name = channel_names[c]
        ax[c].set_xlabel('{0} intensity'.format(channel_name), fontsize = 18)
        ax[c].set_ylabel('# points')
        
        for group in px_groups:
            
            if log_scale:
                ax[c].set_yscale('log')
                
            if colors_dict == None:
                ax[c].hist(px_vals[group][channel_name], bins, label = group, alpha = 0.4)
            else:
                ax[c].hist(px_vals[group][channel_name], bins, label = group, alpha = 0.4, color = colors_dict[group])
            
            
        
        
        if do_bin:
            print('Finding threshold for binarizing')
            if bin_group == None:
                print('Please specify group for binarization')
            else:
                if thresh_scale == None:
                    print('Please specify group for binarization')
                else:
                    all_points = px_vals[bin_group][channel_name]
                    center = np.median(all_points)
                    neg_points = all_points[all_points < center]
                    pos_points = neg_points + center
                    all_points = np.concatenate([neg_points, pos_points])
                    thresh[channel_name] = center + np.std(all_points)*thresh_scale
                    y = np.linspace(ax[c].get_ylim()[0], ax[c].get_ylim()[1], 20)
                    x = np.ones(20)*thresh[channel_name]
                    ax[c].plot(x, y, linestyle = '--', color = 'r', label = 'Binarization threshold')
        
        ax[c].legend(fontsize = 18)

    if not title == None:
        ax[0].set_title(title, fontsize = 18)

    if save:
        if save_loc == None:
            print('Please provide save location')
        else:
            if save_file == None:
                print('Please provide save file name')
            else:
                plt.savefig('{0}{1}{2}'.format(save_loc, sep, save_file))
                
                
    return thresh
