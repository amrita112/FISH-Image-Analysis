"""
Created on Wed 1 Apr 2020

@author: Amrita Singh
"""
import time
import pickle as pkl
from os.path import sep
from PIL import Image
from utils import diff_gauss
from utils import find_threshold
from tifffile import imwrite
from segmentation import get_masks
import numpy as np

def get_bootstrap(data_path, metadata_file, genes = None, n_sample = 1000):

    # Load metadata
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)
    raw_image_path = metadata['raw_image_path']
    filt_image_path = metadata['filt_image_path']
    n_planes = metadata['n_planes']
    h = metadata['h']
    w = metadata['w']
    base_filename = metadata['base_filename']
    channel_names = metadata['channel_names']
    sigma_small = metadata['signal_sigma_small']
    sigma_large = metadata['signal_sigma_large']
    thresh_scale = metadata['signal_thresh_scale']
    if genes == None:
        genes = metadata['genes']

    # Load background pixels
    bg_pixels_no_lipo_file = metadata['bg_pixels_no_lipo_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_pixels_no_lipo_file), 'rb') as f:
        bg_pixels_no_lipo = pkl.load(f)

    cell_pixels_no_lipo_file = metadata['cell_pixels_no_lipo_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_pixels_no_lipo_file), 'rb') as f:
        cell_pixels_no_lipo = pkl.load(f)

    # Allocate variables to store bootstrap values
    bs_raw = {gene: {cell: {} for cell in cell_pixels_no_lipo.keys()} for gene in genes}
    bs_filt = {gene:{cell: {} for cell in cell_pixels_no_lipo.keys()} for gene in genes}

    t0 = time.time()
    for gene in genes:

        print('Loading {0} images: {1} seconds'.format(gene, np.round(time.time() - t0)))

        for plane in range(n_planes):
            img = Image.open('{0}{1}{2}{3}{4}.tif'.format(raw_image_path, sep, base_filename, str(plane + 1).zfill(2),
                                                            channel_names[gene]))
            im_raw = np.array(img)

            img = Image.open('{0}{1}{2}_plane{3}_filt_sigma={4},{5}_thresh={6}.tif'.format(filt_image_path, sep, gene,
                                                                            str(plane + 1).zfill(2),
                                                                            sigma_small, sigma_large,
                                                                            thresh_scale))
            im_filt = np.array(img)
            # Get bootstrap values in raw and filtered images
            cells = metadata['cells_in_plane'][plane]
            for cell in cells:

                bg_px = bg_pixels_no_lipo[cell][plane]
                cell_px = cell_pixels_no_lipo[cell][plane]
                n_verts_bg = bg_px.shape[0]
                n_verts_cell = cell_px.shape[0]
                bs_raw[gene][cell][plane] = 0
                bs_filt[gene][cell][plane] = 0
                if np.logical_or(n_verts_bg == 0, n_verts_cell == 0):
                    print('Cell {0} plane {1}: all pixels excluded from cell or background'.format(cell, plane))
                else:
                    samples = np.random.choice(np.linspace(0, n_verts_bg - 1, n_verts_bg).astype(int),
                                           (n_verts_cell, n_sample), replace = True)

                    for s in range(n_sample):
                        x = bg_px[samples[:, s], 0]
                        y = bg_px[samples[:, s], 1]

                        bs_raw[gene][cell][plane] += np.mean(im_raw[x, y])
                        bs_filt[gene][cell][plane] += np.mean(im_filt[x, y])

                    bs_raw[gene][cell][plane] /= n_sample
                    bs_filt[gene][cell][plane] /= n_sample

    # Save bootstrap values
    bs_raw_file = metadata['bs_raw_file']
    with open('{0}{1}{2}'.format(data_path, sep, bs_raw_file), 'wb') as f:
        pkl.dump(bs_raw, f)

    bs_filt_file = metadata['bs_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, bs_filt_file), 'wb') as f:
        pkl.dump(bs_filt, f)

    return bs_raw, bs_filt

def get_p_values(data_path, metadata_file,  genes = None,):

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

    # Load bootstrap values
    bs_raw_file = metadata['bs_raw_file']
    with open('{0}{1}{2}'.format(data_path, sep, bs_raw_file), 'rb') as f:
        bs_raw = pkl.load(f)

    bs_filt_file = metadata['bs_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, bs_filt_file), 'rb') as f:
        bs_filt = pkl.load(f)

    # Load cell data
    cell_data_file = metadata['cell_data_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_data_file), 'rb') as f:
        cell_data = pkl.load(f)

    cells = list(cell_data.keys())

    pv_raw = {gene: {cell: {} for cell in cells} for gene in genes}
    pv_filt = {gene:{cell: {} for cell in cells} for gene in genes}

    for gene in genes:

        sig = signal_raw[gene]
        sig = np.array([np.mean([sig[cell][plane] for plane in sig[cell].keys()]) for cell in sig.keys()])
        bg = bg_raw[gene]
        bg = np.array([np.mean([bg[cell][plane] for plane in bg[cell].keys()]) for cell in bg.keys()])

        for cell in cells:

            cell_dict = cell_data[cell]
            planes = list(bs_raw[gene][cell].keys())

            bs_means = bs_raw[gene][cell][planes[0]]
            for plane in planes[1:]:
                bs_means = np.concatenate([bs_means, bs_raw[gene][cell][plane]])

            pv_raw[gene][cell] = np.sum(bs_means > sig[cell - 1] - bg[cell - 1])/len(bs_means)

        sig = signal_filt[gene]
        sig = np.array([np.mean([sig[cell][plane] for plane in sig[cell].keys()]) for cell in sig.keys()])
        bg = bg_filt[gene]
        bg = np.array([np.mean([bg[cell][plane] for plane in bg[cell].keys()]) for cell in bg.keys()])

        for cell in cells:

            cell_dict = cell_data[cell]
            planes = cell_dict['z_planes']

            bs_means = bs_filt[gene][cell][planes[0]]
            for plane in planes[1:]:
                bs_means = np.concatenate([bs_means, bs_raw[gene][cell][plane]])

            pv_filt[gene][cell] = np.sum(bs_means > sig[cell - 1] - bg[cell - 1])/len(bs_means)

    # Save p values
    pv_raw_file = metadata['pv_raw_file']
    with open('{0}{1}{2}'.format(data_path, sep, pv_raw_file), 'wb') as f:
        pkl.dump(pv_raw, f)

    pv_filt_file = metadata['pv_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, pv_filt_file), 'wb') as f:
        pkl.dump(pv_filt, f)

    return pv_raw, pv_filt
