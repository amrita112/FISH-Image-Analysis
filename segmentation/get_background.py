"""
Created on Tue 31 Mar 2020

@author: Amrita Singh
"""
import pickle as pkl
import numpy as np
import time
from os.path import sep

def get_bg_pixels(data_path, metadata_file, local_region_um = 50, min_dist_um = 10):

    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)

    try:
        bg_pixels_file = metadata['bg_pixels_file']
        with open('{0}{1}{2}'.format(data_path, sep, bg_pixels_file), 'rb') as f:
            bg_pixels = pkl.load(f)
        print('Background pixels found in {0}{1}{2}'.format(data_path, sep, bg_pixels_file))

        centers_file = metadata['centers_file']
        with open('{0}{1}{2}'.format(data_path, sep, centers_file), 'rb') as f:
            centers = pkl.load(f)

    except IOError:

        cell_data_file = metadata['cell_data_file']
        with open('{0}{1}{2}'.format(data_path, sep, cell_data_file), 'rb') as f:
            cell_data = pkl.load(f)

        cell_pixels_file = metadata['cell_pixels_file']
        with open('{0}{1}{2}'.format(data_path, sep, cell_pixels_file), 'rb') as f:
            cell_pixels = pkl.load(f)

        cells = list(cell_data.keys())

        print('Finding pixels in local background of cells')
        w = metadata['w']
        h = metadata['h']
        n = metadata['n_planes']
        um_per_px = metadata['um_per_px']

        # Find bg_all: to exclude all cell pixels from local background
        bg_all = np.ones([n, h, w])
        n_vert_all = {}
        centers = {}

        for cell in cells:

            cell_dict = cell_data[cell]
            z_planes = cell_dict['z_planes']
            n_vert_all[cell] = {}
            centers[cell] = {}

            for plane in z_planes:

                x = cell_pixels[cell][plane][0, :]
                y = cell_pixels[cell][plane][1, :]

                centers[cell][plane] = np.array([np.mean(x), np.mean(y)]).astype(int)
                bg_all[int(plane), x, y] = np.zeros([1, len(x)])
                n_vert_all[cell][plane] = len(x)

        local_px = int(local_region_um/um_per_px/2)
        min_dist_px = min_dist_um/um_per_px

        bg_pixels = {}

        for cell in cells:
            cell_dict = cell_data[cell]
            z_planes = cell_dict['z_planes']
            bg_pixels[cell] = {}

            p = 0
            for plane in z_planes:

                center = centers[cell][plane]
                n_verts = 4*n_vert_all[cell][plane]

                x1 = np.max([0, center[0] - local_px])
                x2 = np.min([h, center[0] + local_px])
                y1 = np.max([0, center[1] - local_px])
                y2 = np.min([w, center[1] + local_px])

                center_local = [local_px, local_px]

                local_bg = bg_all[int(plane), x1:x2, y1:y2].astype(bool)
                verts = np.array(np.where(local_bg))
                dist_from_center = np.linalg.norm(verts - np.reshape(center_local, [2, 1]), axis = 0)
                verts = verts[:, dist_from_center > min_dist_px]
                dist_from_center = dist_from_center[dist_from_center > min_dist_px]

                order = np.argsort(dist_from_center)
                n_verts = np.min([n_verts, len(order)])


                bg_coords = np.zeros([n_verts, 2])
                bg_coords[:, 0] = verts[0][order[0:n_verts]] + x1
                bg_coords[:, 1] = verts[1][order[0:n_verts]] + y1
                bg_pixels[cell][plane] = bg_coords.astype(int)

        with open('{0}{1}{2}'.format(data_path, sep, bg_pixels_file), 'wb') as f:
            pkl.dump(bg_pixels, f)

        centers_file = metadata['centers_file']
        with open('{0}{1}{2}'.format(data_path, sep, centers_file), 'wb') as f:
            pkl.dump(centers, f)

    return bg_pixels, centers

def get_all_bg_pixels(data_path, metadata_file):

    # Load metadata
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)

    try:
        all_bg_pixels_file = metadata['all_bg_pixels_file']
        with open('{0}{1}{2}'.format(data_path, sep, all_bg_pixels_file), 'rb') as f:
            all_bg_pixels = pkl.load(f)

        print('Background pixels found in {0}{1}{2}'.format(data_path, sep, all_bg_pixels_file))

    except IOError:
        # Load cell pixels
        bg_pixels_file = metadata['bg_pixels_file']
        with open('{0}{1}{2}'.format(data_path, sep, bg_pixels_file), 'rb') as f:
            bg_pixels = pkl.load(f)

        # Load cell data
        cell_data_file = metadata['cell_data_file']
        with open('{0}{1}{2}'.format(data_path, sep, cell_data_file), 'rb') as f:
            cell_data = pkl.load(f)

        cells = bg_pixels.keys()
        all_bg_pixels = np.zeros([1, 3])

        for cell in cells:
            planes = cell_data[cell]['z_planes']
            for plane in planes:
                xvals = np.reshape(bg_pixels[cell][plane][0], [-1, 1])
                yvals = np.reshape(bg_pixels[cell][plane][1], [-1, 1])
                zvals = np.ones([len(xvals), 1])*int(plane)
                coords = np.concatenate((zvals, xvals, yvals), axis = 1)
                all_bg_pixels = np.concatenate((all_bg_pixels, coords), axis = 0)

        all_bg_pixels = all_bg_pixels[1:, :].astype(int)
        with open('{0}{1}{2}'.format(data_path, sep, all_bg_pixels_file), 'wb') as f:
            pkl.dump(all_bg_pixels, f)

    return all_bg_pixels
