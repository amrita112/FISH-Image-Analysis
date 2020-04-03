"""
Created on Wed 1 Apr 2020

@author: Amrita Singh
"""
import pickle as pkl
from utils import pixels_in_roi
from os.path import sep

def get_rois(data_path, metadata_file):
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)
    lipo_rois_file = metadata['lipo_rois_file']
    with open('{0}{1}{2}'.format(data_path, sep, lipo_rois_file), 'rb') as f:
        l_rois = pkl.load(f)
    print('{0} rois'.format(len(l_rois)))

    return l_rois

def get_cell_lipo_pixels(data_path, metadata_file):
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)
    lipo_pixels_cells_file = metadata['lipo_pixels_cells_file']
    try:
        with open('{0}{1}{2}'.format(data_path, sep, lipo_pixels_cells_file), 'rb') as f:
            data_dict = pkl.load(f)
        all_lipo_pixels = data_dict['all_lipo_pixels']
        lipo_pixels_roi = data_dict['lipo_pixels_roi']
        print('Data loaded')

    except IOError:
        print('No saved data found, finding pixels in lipofuscin rois')

        h = metadata['h']
        w = metadata['w']
        n_planes = metadata['n_planes']

        lipo_rois_file = metadata['lipo_rois_file']
        with open('{0}{1}{2}'.format(data_path, sep, lipo_rois_file), 'rb') as f:
            l_rois = pkl.load(f)

        px = pixels_in_roi.pixels_in_roi(h, w, n_planes, l_rois)
        all_lipo_pixels = px['all_pixels']
        lipo_pixels_roi = px['pixels_roi']
        with open('{0}{1}{2}'.format(data_path, sep, lipo_pixels_cells_file), 'wb') as f:
            pkl.dump({'all_lipo_pixels': all_lipo_pixels, 'lipo_pixels_roi': lipo_pixels_roi}, f)

    return all_lipo_pixels, lipo_pixels_roi
