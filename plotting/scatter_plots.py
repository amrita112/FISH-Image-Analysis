"""
Created on Wed 1 Apr 2020

@author: Amrita Singh
"""
import time
import pickle as pkl
from os.path import sep
import numpy as np
import matplotlib.pyplot as plt

def signal_vs_depth(data_path, metadata_file, genes = None, save = False, title = None,):

    # Load metadata
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)
    if genes == None:
        genes = metadata['genes']

    # Load signal
    signal_raw_file = metadata['signal_raw_file']
    with open('{0}{1}{2}'.format(data_path, sep, signal_raw_file), 'rb') as f:
        signal_raw = pkl.load(f)

    signal_filt_file = metadata['signal_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, signal_filt_file), 'rb') as f:
        signal_filt = pkl.load(f)

    # Load background
    bg_raw_file = metadata['bg_raw_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_raw_file), 'rb') as f:
        bg_raw = pkl.load(f)

    bg_filt_file = metadata['bg_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_filt_file), 'rb') as f:
        bg_filt = pkl.load(f)

    # Load depths
    depths_file = metadata['depths_file']
    with open('{0}{1}{2}'.format(data_path, sep, depths_file), 'rb') as f:
        d = pkl.load(f)

    # Plot intensity vs depth for all channels
    fig, ax = plt.subplots(nrows = 2, ncols = len(genes), figsize = (20, 15))

    for g in range(len(genes)):
        gene = genes[g]

        sig = signal_raw[gene]
        sig = np.array([np.mean([sig[cell][plane] for plane in sig[cell].keys()]) for cell in sig.keys()])
        bg = bg_raw[gene]
        bg = np.array([np.mean([bg[cell][plane] for plane in bg[cell].keys()]) for cell in bg.keys()])
        ax[0, g].scatter(d, sig - bg, color = 'k', alpha = 0.5, marker = '.')
        ax[0, g].set_xlabel('Depth from pia (um)')
        ax[0, g].set_ylabel('Signal - background in raw images')
        ax[0, g].set_title(gene)

        sig = signal_filt[gene]
        sig = np.array([np.mean([sig[cell][plane] for plane in sig[cell].keys()]) for cell in sig.keys()])
        bg = bg_filt[gene]
        bg = np.array([np.mean([bg[cell][plane] for plane in bg[cell].keys()]) for cell in bg.keys()])
        ax[1, g].scatter(d, sig - bg, color = 'k', alpha = 0.5, marker = '.')
        ax[1, g].set_xlabel('Depth from pia (um)')
        ax[1, g].set_ylabel('Signal - background in filtered images')
        ax[1, g].set_title(gene)
