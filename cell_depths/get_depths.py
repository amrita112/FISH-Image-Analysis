"""
Created on Wed 1 Apr 2020

@author: Amrita Singh
"""
import pickle as pkl
from os.path import sep
from plotting import napari_rendering
import numpy as np

def draw_surface(data_path, metadata_file, viewer):

    surface_layer = napari_rendering.draw_surface(data_path, metadata_file, viewer)
    return surface_layer

def approx_altitude(d1, d2, cell_center):
    """d1 and d2 are two points on the surface closest to the cell center"""

    dx = d2[0] - d1[0]
    dy = d2[1] - d1[1]
    d1_cell_x = cell_center[0] - d1[0]
    d1_cell_y = cell_center[1] - d1[1]
    t = (dx * d1_cell_x + dy * d1_cell_y) / (dx * dx + dy * dy)
    px = d1[0] + t * dx # Projection of cell center onto line joining d1 and d2
    py = d1[1] + t * dy
    depth = np.linalg.norm([px - cell_center[0], py - cell_center[1]])
    return depth

def get_depths(data_path, metadata_file, surface_layer = None):

    # Load metadata
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)

    # Check if depths already exist
    try:
        depths_file = metadata['depths_file']
        with open('{0}{1}{2}'.format(data_path, sep, depths_file), 'rb') as f:
            depths = pkl.load(f)
        print('Depths loaded')

    except:
        print('Could not load depths. Draw surface to calculate depths')
        # Load cell masks
        cell_data_file = metadata['cell_data_file']
        with open('{0}{1}{2}'.format(data_path, sep, cell_data_file), 'rb') as f:
            cell_data = pkl.load(f)

        # Load cell centers
        centers_file = metadata['centers_file']
        with open('{0}{1}{2}'.format(data_path, sep, centers_file), 'rb') as f:
            centers = pkl.load(f)

        depths_file = metadata['depths_file']

        if surface_layer == None:
            print('Draw surface layer in Napari')
            return
        
        surface = surface_layer.data[-1][:, 1:]
        cells = list(cell_data.keys())
        depths = np.zeros(len(cells))
        um_per_px = metadata['um_per_px']

        for cell in cells:
            pls = centers[cell].keys()
            cx = []
            cy = []
            for pl in pls:
                cx = np.append(cx, centers[cell][pl][0])
                cy = np.append(cy, centers[cell][pl][1])
            c = [np.mean(cx), np.mean(cy)]
            distances = np.linalg.norm(c - surface, axis = 1)
            order = np.argsort(distances)
            [p1, p2] = order[0:2]
            d1 = surface[p1]
            d2 = surface[p2]
            depths[cell - 1] = approx_altitude(d1, d2, c)
        depths = depths*um_per_px

        depths_file = metadata['depths_file']
        with open('{0}{1}{2}'.format(data_path, sep, depths_file), 'wb') as f:
            pkl.dump(depths, f)

    return depths
