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

def filter(data_path, metadata_file, genes = None, sigma_small = 5, sigma_large = 10, thresh_scale = 1,
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

    cell_pixels_file = metadata['cell_pixels_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_pixels_file), 'rb') as f:
        cell_pixels = pkl.load(f)

    # Load background pixels
    bg_pixels_file = metadata['bg_pixels_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_pixels_file), 'rb') as f:
        bg_pixels = pkl.load(f)

    # Load annotated lipofuscin pixels
    lipo_pixels_cells_file = metadata['lipo_pixels_cells_file']
    with open('{0}{1}{2}'.format(data_path, sep, lipo_pixels_cells_file), 'rb') as f:
        data_dict = pkl.load(f)
        all_lipo_pixels = data_dict['all_lipo_pixels']

    # Allocate variables to store pixel values in filtered images
    cell_filt_values = {gene: {cell: {} for cell in cell_pixels.keys()} for gene in genes}
    bg_filt_values = {gene:{cell: {} for cell in cell_pixels.keys()} for gene in genes}
    annotated_lipo_filt_values = {gene: [] for gene in genes}

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

            # Get pixel values for cells and background in filtered images
            cells_in_plane = get_masks.get_cells_in_plane(data_path, metadata_file, plane)
            for cell in cells_in_plane:

                cell_px = cell_pixels[cell][plane]
                cell_filt_values[gene][cell][plane] = im_thresh[cell_px[0, :], cell_px[1, :]]

                bg_px = bg_pixels[cell][plane]
                bg_filt_values[gene][cell][plane] = im_thresh[bg_px[:, 0], bg_px[:, 1]]

            # Get pixel values for annotated lipofuscin pixels in filtered images
            lipo_px = all_lipo_pixels[all_lipo_pixels[:, 0] == plane, 1:]
            annotated_lipo_filt_values[gene] = np.append(annotated_lipo_filt_values[gene],
                                                            im_thresh[lipo_px[:, 0], lipo_px[:, 1]])


    # Save cell pixel values
    cell_pixels_filt_file = metadata['cell_pixels_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_pixels_filt_file), 'wb') as f:
        pkl.dump(cell_filt_values, f)

    # Save background pixel values
    bg_pixels_filt_file = metadata['bg_pixels_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_pixels_filt_file), 'wb') as f:
        pkl.dump(bg_filt_values, f)

    # Load annotated lipofuscin pixels
    lipo_pixels_cells_filt_file = metadata['lipo_pixels_cells_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, lipo_pixels_cells_filt_file), 'wb') as f:
        pkl.dump(annotated_lipo_filt_values, f)

    metadata['lipo_sigma_small'] = sigma_small
    metadata['lipo_sigma_large'] = sigma_large
    metadata['lipo_thresh_scale'] = thresh_scale
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'wb') as f:
        pkl.dump(metadata, f)

    return cell_filt_values, bg_filt_values, annotated_lipo_filt_values
