# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 15:36:38 2020

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


def scatter_channels(img_dict, pixels_dict, max_per_group = 10000, colors_dict = None,
                     make_3D = False, alpha_dict = None,
                     title = None, save = False, save_loc = None, save_file = None, log_scale = False):


    # %%
    channel_names = list(img_dict.keys())
    n_channels = len(channel_names)
    channels = range(n_channels)

    if make_3D:
        combos = list(combinations(channels, 3))
    else:
        combos = list(combinations(channels, 2))

    n_combos = len(combos)

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
            if type(colors_dict[group]) is np.ndarray:
                if colors_dict[group].shape[0] > max_per_group:
                    colors_dict[group] = colors_dict[group][idx == 0]

        px_vals[group] = {}

        for channel in channel_names:
            px_vals[group][channel] = img_dict[channel][(px[:, 0], px[:, 1], px[:, 2])]

    # %%
    fig, ax = plt.subplots(nrows = n_combos if n_combos > 1 else 2, ncols = 1, figsize = (20, 20))
    for combo in range(n_combos):

        if make_3D:
            c3 = channel_names[combos[combo][2]]
            ax[combo] = plt.axes(projection='3d')
            ax[combo].set_zlabel(c3, fontsize = 18)

        c1 = channel_names[combos[combo][0]]
        c2 = channel_names[combos[combo][1]]

        ax[combo].set_xlabel(c1, fontsize = 18)
        ax[combo].set_ylabel(c2, fontsize = 18)



        for group in px_groups:

            if log_scale:
                ax[combo].set_xscale('log')
                ax[combo].set_yscale('log')
                if make_3D:
                    ax[combo].set_zscale('log')

            if colors_dict == None:
                if make_3D:
                    p = ax[combo].scatter3D(px_vals[group][c1], px_vals[group][c2], px_vals[group][c3], label = group,
                      marker = '.', alpha = alpha_dict[group])
                    fig.colorbar(p)
                else:
                    ax[combo].scatter(px_vals[group][c1], px_vals[group][c2], label = group, marker = '.', alpha = alpha_dict[group])
            else:
                if make_3D:
                    ax[combo].scatter3D(px_vals[group][c1], px_vals[group][c2], px_vals[group][c3], label = group,
                      marker = '.', alpha = alpha_dict[group], c = colors_dict[group])
                else:
                    ax[combo].scatter(px_vals[group][c1], px_vals[group][c2], label = group, c = colors_dict[group], marker = '.', alpha = alpha_dict[group])


        ax[combo].legend(fontsize = 18)


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
