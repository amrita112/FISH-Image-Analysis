"""
Created on Tue 31 Mar 2020

@author: Amrita Singh
"""
import pickle as pkl
import numpy as np
import time
import matplotlib.path as mpltpath
from os.path import sep

def get_masks(data_path, metadata_file, ):
    """Returns  cell_data as a dictionary and n_cells """
    #print('Hello world')
    cell_data = None
    n_cells = None
    try:
        with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
            metadata = pkl.load(f)
        cell_data_file = metadata['cell_data_file']
        with open('{0}{1}{2}'.format(data_path, sep, cell_data_file), 'rb') as f:
            cell_data = pkl.load(f)
        indices = list(cell_data.keys())
        if not np.max(indices) == len(indices):
            print('Re-numbering cells to be consecutive')
            cell_data_temp = {}
            for i in range(len(indices)):
                cell_data_temp[i + 1] = cell_data[indices[i]]
                cell_data_temp[i + 1]['cell_id'] = i + 1
            cell_data = cell_data_temp
            with open('{0}{1}{2}'.format(data_path, sep, data_file), 'wb') as f:
                pkl.dump(cell_data, f)
            n_cells = i + 1
        else:
            n_cells = len(indices)
        print('{0} cell masks found'.format(n_cells))

    except:
        print('No data found')
    return cell_data, n_cells

def get_cell_pixels(data_path, metadata_file):

    t0 = time.time()

    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)

    try:
        cell_pixels_file = metadata['cell_pixels_file']
        with open('{0}{1}{2}'.format(data_path, sep, cell_pixels_file), 'rb') as f:
            cell_pixels = pkl.load(f)

        print('Cell pixels found in {0}{1}{2}'.format(data_path, sep, cell_pixels_file))

    except IOError:
        print('Finding pixels in cells')
        w = metadata['w']
        h = metadata['h']

        cell_data_file = metadata['cell_data_file']
        with open('{0}{1}{2}'.format(data_path, sep, cell_data_file), 'rb') as f:
            cell_data = pkl.load(f)

        xv = range(w)
        yv = range(h)
        coord_array = np.array(np.meshgrid(xv, yv))
        points = np.zeros([h*w, 2])
        p = 0
        for i in range(h):
            for j in range(w):
                points[p, 1] = coord_array[0, i, j]
                points[p, 0] = coord_array[1, i, j]
                p += 1

        no_cells = len(cell_data.keys())
        cell_pixels = {}
        for cell in range(no_cells):
            if np.mod(cell, 2) == 0:
                print('Cell {0}: {1} seconds'.format(cell, np.round(time.time() - t0)))
            cell_no = cell + 1
            cell_pixels[cell_no] = {}
            cell_dict = cell_data[cell_no]
            masks = cell_dict['masks']
            z_planes = cell_dict['z_planes']
            for plane in z_planes:

                vertices = masks[plane]
                path = mpltpath.Path(vertices)
                mask = path.contains_points(points)
                mask = np.reshape(mask, [h, w])
                cell_pixels[cell_no][plane] = np.array(np.where(mask))

        with open('{0}{1}{2}'.format(data_path, sep, cell_pixels_file), 'wb') as f:
            pkl.dump(cell_pixels, f)

    return cell_pixels


def get_all_cell_pixels(data_path, metadata_file):

    # Load metadata
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)

    try:
        all_cell_pixels_file = metadata['all_cell_pixels_file']
        with open('{0}{1}{2}'.format(data_path, sep, all_cell_pixels_file), 'rb') as f:
            all_cell_pixels = pkl.load(f)

        print('Cell pixels found in {0}{1}{2}'.format(data_path, sep, all_cell_pixels_file))

    except IOError:
        # Load cell pixels
        cell_pixels_file = metadata['cell_pixels_file']
        with open('{0}{1}{2}'.format(data_path, sep, cell_pixels_file), 'rb') as f:
            cell_pixels = pkl.load(f)

        # Load cell data
        cell_data_file = metadata['cell_data_file']
        with open('{0}{1}{2}'.format(data_path, sep, cell_data_file), 'rb') as f:
            cell_data = pkl.load(f)

        cells = cell_pixels.keys()
        all_cell_pixels = np.zeros([1, 3])

        for cell in cells:
            planes = cell_data[cell]['z_planes']
            for plane in planes:
                xvals = np.reshape(cell_pixels[cell][plane][0], [-1, 1])
                yvals = np.reshape(cell_pixels[cell][plane][1], [-1, 1])
                zvals = np.ones([len(xvals), 1])*plane
                coords = np.concatenate((zvals, xvals, yvals), axis = 1)
                all_cell_pixels = np.concatenate((all_cell_pixels, coords), axis = 0)

        all_cell_pixels = all_cell_pixels[1:, :].astype(int)
        with open('{0}{1}{2}'.format(data_path, sep, all_cell_pixels_file), 'wb') as f:
            pkl.dump(all_cell_pixels, f)

    return all_cell_pixels

def get_cells_in_plane(data_path, metadata_file, query_plane):
    """ Return cells with a mask in the plane """
    try:
        with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
            metadata = pkl.load(f)
        cells_in_plane = metadata['cells_in_plane']
    except:

        cells_in_plane = {}

        cell_data_file = metadata['cell_data_file']
        with open('{0}{1}{2}'.format(data_path, sep, cell_data_file), 'rb') as f:
            cell_data = pkl.load(f)

        cells = list(cell_data.keys())
        for cell in cells:
            planes = cell_data[cell]['z_planes']
            for plane in planes:
                if plane in cells_in_plane.keys():
                    cells_in_plane[plane] = np.append(cells_in_plane[plane], cell)
                else:
                    cells_in_plane[plane] = [cell]

        metadata['cells_in_plane'] = cells_in_plane
        with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'wb') as f:
            pkl.dump(metadata, f)

    return cells_in_plane[query_plane]
