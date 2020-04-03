"""
Created on Tue 31 Mar 2020

@author: Amrita Singh
"""

import napari
from os.path import sep
import pickle as pkl
import numpy as np
import time
from PIL import Image

def create_viewer():

    viewer = napari.Viewer()
    return viewer

def images(data_path, metadata_file, viewer, planes = None, genes = None, colors = ['red', 'green', 'blue', 'magenta'],
            image_type = 'raw'):
    """ Render raw images in napari viewer. Image type = 'raw' (default) or 'filt'.
    Default colors: ['red', 'green', 'blue', 'magenta']"""
    t0 = time.time()
    # Load metadata
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)
    if image_type == 'raw':
        image_path = metadata['raw_image_path']
    else:
        image_path = metadata['filt_image_path']
    if genes == None:
        genes = metadata['genes']
    if planes == None:
        planes = list(range(metadata['n_planes']))
    n_planes = len(planes)
    channel_names = metadata['channel_names']
    base_filename = metadata['base_filename']
    h = metadata['h']
    w = metadata['w']
    if len(genes) > 4:
        print('More than 4 channels - please specify colors to use for display')
        return

    im_array = np.zeros([n_planes, h, w])
    for g in range(len(genes)):
        gene = genes[g]
        print('Loading {0} images: {1} seconds'.format(gene, np.round(time.time() - t0)))

        for plane in planes:
            print('     Plane {0}, {1} seconds'.format(plane + 1, np.round(time.time() - t0)))
            img = Image.open('{0}{1}{2}{3}{4}.tif'.format(image_path, sep, base_filename, str(plane + 1).zfill(2),
                                                            channel_names[gene]))
            im_array[plane, :, :] = np.array(img)

        viewer.add_image(im_array, name = '{0} {1}'.format(gene, image_type),
                            colormap = colors[g], blending = 'additive')

def filt_images(data_path, metadata_file, viewer, planes = None, genes = None, colors = ['red', 'green', 'blue', 'magenta'],
            image_type = 'filt'):
    """ Render raw images in napari viewer. Image type = 'raw' (default) or 'filt'.
    Default colors: ['red', 'green', 'blue', 'magenta']"""
    t0 = time.time()
    # Load metadata
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)
    if image_type == 'filt':
        image_path = metadata['filt_image_path']
    if genes == None:
        genes = metadata['genes']
    if planes == None:
        planes = list(range(metadata['n_planes']))
    n_planes = len(planes)
    channel_names = metadata['channel_names']
    base_filename = metadata['base_filename']
    sigma_small = metadata['lipo_sigma_small']
    sigma_large = metadata['lipo_sigma_large']
    thresh_scale = metadata['lipo_thresh_scale']
    h = metadata['h']
    w = metadata['w']
    if len(genes) > 4:
        print('More than 4 channels - please specify colors to use for display')
        return

    im_array = np.zeros([n_planes, h, w])
    for g in range(len(genes)):
        gene = genes[g]
        print('Loading {0} images: {1} seconds'.format(gene, np.round(time.time() - t0)))

        for plane in planes:
            print('     Plane {0}, {1} seconds'.format(plane + 1, np.round(time.time() - t0)))
            img = Image.open('{0}{1}{2}_plane{3}_filt_sigma={4},{5}_thresh={6}.tif'.format(image_path, sep, gene,
                                                                                        str(plane + 1).zfill(2),
                                                                                        sigma_small, sigma_large,
                                                                                        thresh_scale))
            im_array[plane, :, :] = np.array(img)

        viewer.add_image(im_array, name = '{0} {1}'.format(gene, image_type),
                            colormap = colors[g], blending = 'additive')

def cell_masks(data_path, metadata_file, viewer, planes = None):

    # Load metadata
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)
    cells_in_plane = metadata['cells_in_plane']

    cell_data_file = metadata['cell_data_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_data_file), 'rb') as f:
        cell_data = pkl.load(f)

    if planes == None:
        planes = range(metadata['n_planes'])

    mask_layer = viewer.add_shapes(data = None, name = 'Cell masks', edge_color = 'white')

    cells = np.concatenate([cells_in_plane[plane] for plane in planes])
    for cell in cells:
         if np.mod(cell, 10) == 0:
             print('Cell {0}'.format(cell))
         cell_planes = cell_data[cell]['z_planes']
         cell_planes = [plane for plane in cell_planes if plane in planes]
         for plane in cell_planes:
             mask = cell_data[cell]['masks'][plane]
             mask = np.concatenate((np.ones([mask.shape[0], 1])*plane, mask), axis = 1)
             mask = np.append(mask, mask[np.newaxis, 0, :], axis = 0)
             mask_layer.add(mask, shape_type = 'path', opacity = 1, edge_width = 3)

def lipo_rois(data_path, metadata_file, viewer, planes = None):

    # Load metadata
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)

    # Load annotated lipofuscin rois
    lipo_rois_file = metadata['lipo_rois_file']
    with open('{0}{1}{2}'.format(data_path, sep, lipo_rois_file), 'rb') as f:
        l_rois = pkl.load(f)

    if planes == None:
        planes = range(metadata['n_planes'])

    lipo_roi_layer = viewer.add_shapes(data = None, name = 'Lipofuscin rois', edge_color = 'white')
    flag = 0

    for roi in l_rois:
        if roi[0, 0] in planes:
             flag = 1
             roi = np.append(roi, roi[np.newaxis, 0, :], axis = 0)
             lipo_roi_layer.add(roi, shape_type = 'path', opacity = 1, edge_width = 3)
    if flag == 0:
        print('No lipofuscin rois in the specified planes')



def no_lipo_masks(data_path, metadata_file, viewer, planes = None):
    # Load metadata
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)
    cells_in_plane = metadata['cells_in_plane']

    # Load cell masks
    cell_data_file = metadata['cell_data_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_data_file), 'rb') as f:
        cell_data = pkl.load(f)

    # Load cell masks excluding Lipofuscin
    cell_pixels_no_lipo_file = metadata['cell_pixels_no_lipo_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_pixels_no_lipo_file), 'rb') as f:
        cell_pixels_no_lipo = pkl.load(f)

    if planes == None:
        planes = range(metadata['n_planes'])

    im_array = np.zeros([len(planes), metadata['h'], metadata['w']])

    for plane in planes:
        cells = metadata['cells_in_plane'][plane]
        for cell in cells:
            px = cell_pixels_no_lipo[cell][plane]
            im_array[np.ones(px.shape[1]).astype(int)*plane, px[0, :], px[1, :]] = np.ones(px.shape[1])

    viewer.add_image(data = im_array, name = 'Cell masks excluding lipo', colormap = 'gray', blending = 'additive')

def draw_surface(data_path, metadata_file, viewer):
    surface_layer = viewer.add_shapes(data = None, name = 'Surface')
    return surface_layer
