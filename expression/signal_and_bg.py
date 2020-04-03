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

def get_signal_and_bg(data_path, metadata_file, genes = None, sigma_small = 1, sigma_large = 2, thresh_scale = 12,
            save_tiffs = False):

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
    if genes == None:
        genes = metadata['genes']

    # Load cell pixels
    all_cell_pixels_file = metadata['all_cell_pixels_file']
    with open('{0}{1}{2}'.format(data_path, sep, all_cell_pixels_file), 'rb') as f:
        all_cell_pixels = pkl.load(f)

    cell_pixels_no_lipo_file = metadata['cell_pixels_no_lipo_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_pixels_no_lipo_file), 'rb') as f:
        cell_pixels_no_lipo = pkl.load(f)

    # Load background pixels
    bg_pixels_no_lipo_file = metadata['bg_pixels_no_lipo_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_pixels_no_lipo_file), 'rb') as f:
        bg_pixels_no_lipo = pkl.load(f)

    # Allocate variables to store pixel values in filtered images
    signal_raw = {gene: {cell: {} for cell in cell_pixels_no_lipo.keys()} for gene in genes}
    signal_filt = {gene:{cell: {} for cell in cell_pixels_no_lipo.keys()} for gene in genes}
    bg_raw = {gene: {cell: {} for cell in cell_pixels_no_lipo.keys()} for gene in genes}
    bg_filt = {gene:{cell: {} for cell in cell_pixels_no_lipo.keys()} for gene in genes}

    t0 = time.time()
    im_thresh = np.zeros([h, w])
    for gene in genes:

        print('Filtering {0} images: {1} seconds'.format(gene, np.round(time.time() - t0)))

        for plane in range(n_planes):
            print('Plane {0}, {1} seconds'.format(plane + 1, np.round(time.time() - t0)))
            img = Image.open('{0}{1}{2}{3}{4}.tif'.format(raw_image_path, sep, base_filename, str(plane + 1).zfill(2),
                                                            channel_names[gene]))
            im_array = np.array(img)
            im_diff = diff_gauss.diff_gauss(sigma_small, sigma_large, im_array)

            px = all_cell_pixels[all_cell_pixels[:, 0] == plane, 1:]
            points = im_diff[px[:, 0], px[:, 1]]
            thresh = find_threshold.find_threshold(points, thresh_scale)

            ix = np.where(im_diff > thresh)
            im_thresh[ix] = im_diff[ix]


            if save_tiffs:
                im_thresh_save = (im_thresh*65536/np.max(im_thresh)).astype(np.int16)
                print('Saving plane {0}: {1} seconds'.format(plane, np.round(time.time() - t0)))
                imwrite('{0}{1}{2}_plane{3}_filt_sigma={4},{5}_thresh={6}.tif'.format(filt_image_path, sep, gene,
                                                                                str(plane + 1).zfill(2),
                                                                                sigma_small, sigma_large,
                                                                                thresh_scale),
                                                                                im_thresh_save)

            # Get mean values for cells and background in raw and filtered images
            cells = metadata['cells_in_plane'][plane]
            for cell in cells:

                cell_px = cell_pixels_no_lipo[cell][plane]
                if cell_px.shape[1] == 0:
                    signal_raw[gene][cell][plane] = 0
                    signal_filt[gene][cell][plane] = 0
                    print('Cell {0} plane {1}: all pixels excluded'.format(cell, plane))
                else:
                    signal_raw[gene][cell][plane] = np.mean(im_array[cell_px[0, :], cell_px[1, :]])
                    signal_filt[gene][cell][plane] = np.mean(im_thresh[cell_px[0, :], cell_px[1, :]])

                bg_px = bg_pixels_no_lipo[cell][plane]
                if bg_px.shape[1] == 0:
                    bg_raw[gene][cell][plane] = 0
                    bg_filt[gene][cell][plane] = 0
                    print('Cell {0} plane {1}: all background pixels excluded'.format(cell, plane))
                else:
                    bg_raw[gene][cell][plane] = np.mean(im_array[bg_px[:, 0], bg_px[:, 1]])
                    bg_filt[gene][cell][plane] = np.mean(im_thresh[bg_px[:, 0], bg_px[:, 1]])


    # Save signal and background values
    signal_raw_file = metadata['signal_raw_file']
    with open('{0}{1}{2}'.format(data_path, sep, signal_raw_file), 'wb') as f:
        pkl.dump(signal_raw, f)

    signal_filt_file = metadata['signal_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, signal_filt_file), 'wb') as f:
        pkl.dump(signal_filt, f)

    # Save background pixel values
    bg_filt_file = metadata['bg_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_filt_file), 'wb') as f:
        pkl.dump(bg_filt, f)

    bg_raw_file = metadata['bg_raw_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_raw_file), 'wb') as f:
        pkl.dump(bg_raw, f)

    metadata['signal_sigma_small'] = sigma_small
    metadata['signal_sigma_large'] = sigma_large
    metadata['signal_thresh_scale'] = thresh_scale
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'wb') as f:
        pkl.dump(metadata, f)

    return signal_raw, signal_filt, bg_raw, bg_filt
