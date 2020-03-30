# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 16:36:35 2020

Segment cells and remove lipofuscin pixels

@author: Amrita S
"""
from IPython import get_ipython
ipython = get_ipython()
ipython.magic("gui qt5") 

import napari
from PIL import Image

import numpy as np
import matplotlib.pyplot as plt
import cv2
import matplotlib.colors as colors
import time

import pickle as pkl
from os.path import sep

import Cell
import pixels_in_roi
import diff_gauss
import scatter_channels
import hist_channels
import m_dist

# %% Load images
t0 = time.time()
folder = 'G:\\Shared drives\\as_share\\HCR\\HCR_10.17\\S2_DAPI_546_647_514_594_2019_10_19__16_28_49'
plane_nos = range(1, 16)

# base_filename2 = 'S1_dapi_546_488_647_s2z'
base_filename = 'S2_DAPI_546_647_514_594_2019_10_19__16_28_49_z'
n = len(plane_nos)
print('Number of planes: {0}'.format(n))

# Create 4D array to store images
img = Image.open('{0}\\{1}{2}c4_ORG.tif'.format(folder, base_filename, str(plane_nos[0]).zfill(2)))

h = img.height
w = img.width
im_array_gad1 = np.zeros([n, h, w])
im_array_vip = np.zeros([n, h, w])
im_array_sst = np.zeros([n, h, w])
im_array_ndnf = np.zeros([n, h, w])
print('Size of image in pixels: {0} X {1} X {2}'.format(n, h, w))

for i in range(n):
    print('Loading image {0}, {1} seconds'.format(i + 1, np.round(time.time() - t0)))
    p = plane_nos[i]
    img_gad1 = Image.open('{0}\{1}{2}c4_ORG.tif'.format(folder, base_filename, str(p).zfill(2)))
    img_vip = Image.open('{0}\{1}{2}c2_ORG.tif'.format(folder, base_filename, str(p).zfill(2)))
    img_sst = Image.open('{0}\{1}{2}c5_ORG.tif'.format(folder, base_filename, str(p).zfill(2)))
    img_ndnf = Image.open('{0}\{1}{2}c3_ORG.tif'.format(folder, base_filename, str(p).zfill(2)))
    try:
        im_array_gad1[i, :, :] = np.array(img_gad1)
        im_array_vip[i, :, :] = np.array(img_vip)
        im_array_sst[i, :, :] = np.array(img_sst)
        im_array_ndnf[i, :, :] = np.array(img_ndnf)
    except:
        print('Plane {0} could not be loaded'.format(p))
        print('Size of plane {0} in pixels: {0} X {1}'.format(img.height, img.width))
        im_array_gad1 = np.delete(im_array_gad1, i, axis = 0)
        im_array_vip = np.delete(im_array_vip, i, axis = 0)
        im_array_sst = np.delete(im_array_sst, i, axis = 0)
        im_array_ndnf = np.delete(im_array_ndnf, i, axis = 0)
        plane_nos.remove(p)
        i -= 1
        n = len(plane_nos)
        continue

del img_gad1
del img_vip
del img_sst
del img_ndnf

# %% Render image in napari gui

viewer = napari.Viewer()
viewer.add_image(im_array_gad1, name = 'Gad1', colormap = 'cyan', blending = 'additive')
viewer.add_image(im_array_ndnf, name = 'Ndnf', colormap = 'magenta', blending = 'additive')
viewer.add_image(im_array_vip, name = 'Vip', colormap = 'yellow', blending = 'additive')
viewer.add_image(im_array_sst, name = 'Sst', colormap = 'green', blending = 'additive')

# %% Load cell rois
# Load masks if they already exist
save_loc = 'G:\\Shared drives\\as_share\\HCR\\HCR_10.17'
save_file= 'S2_data.pkl'
try:
    with open('{0}\{1}'.format(save_loc, save_file), 'rb') as f:
        Cell.cell_data = pkl.load(f)
        indices = list(Cell.cell_data.keys())
        if not np.max(indices) == len(indices):
            print('Re-numbering cells to be consecutive')
            Cell.cell_data_temp = {}
            for i in range(len(indices)):
                Cell.cell_data_temp[i + 1] = Cell.cell_data[indices[i]]
                Cell.cell_data_temp[i + 1]['cell_id'] = i + 1
            Cell.cell_data = Cell.cell_data_temp 
            with open('{0}\{1}'.format(save_loc, save_file), 'wb') as f:
                pkl.dump(Cell.cell_data, f)
            Cell.n_cells = i + 1
        else:
            Cell.n_cells = len(indices)
    print('{0} cells found'.format(Cell.n_cells))
except:
    print('No data found')


# %% Add masks to napari viewer

mask_layer = viewer.add_shapes(name = 'Cell masks')

indices = list(Cell.cell_data.keys())
for cell in indices:
    if np.mod(cell, 10) == 0:
        print('Cell {0}'.format(cell))
    planes = Cell.cell_data[cell]['z_planes']
    for plane in planes:
        mask = Cell.cell_data[cell]['masks'][plane]
        mask = np.concatenate((np.ones([mask.shape[0], 1])*plane, mask), axis = 1)
        mask_layer.add(mask, shape_type = 'polygon', opacity = 0.2, face_color = 'white', edge_color = 'red', edge_width = 3)


# %% Load lipofuscin rois
save_loc = 'G:\\Shared drives\\as_share\\HCR\\HCR_10.17'
save_file = 'HCR_10.17_S2_lipofuscin_rois_in_cells.pkl'

with open('{0}\\{1}'.format(save_loc, save_file), 'rb') as f:
    l_rois = pkl.load(f)    

# %% Add lipofuscin rois to viewer
n_rois = len(l_rois)
print('{0} rois found'.format(n_rois))
viewer.add_shapes(l_rois, name = 'Lipofuscin ROIs',
                   shape_type = 'polygon', opacity = 1, face_color = 'white', 
                   edge_color = 'blue', edge_width = 3)

# %% Save lipofuscin rois from viewer
save_loc = 'G:\\Shared drives\\as_share\\HCR\\HCR_10.17'
save_file = 'HCR_10.17_S2_lipofuscin_rois_in_cells.pkl'

l_rois = viewer.layers['Lipofuscin ROIs'].data

with open('{0}\\{1}'.format(save_loc, save_file), 'wb') as f:
    pkl.dump(l_rois, f)    

# %% Get pixels in cell ROIs
save_loc = 'G:\\Shared drives\\as_share\\HCR\\HCR_10.17'
# save_file = 'S2_cell_pixels.pkl'
save_file = 'S2_mask_vertices.pkl'

try:
    with open('{0}\\{1}'.format(save_loc, save_file), 'rb') as f:
        mask_vertices = pkl.load(f)
        cells = mask_vertices.keys()
        x = 0
        all_cell_pixels = np.zeros([1, 3])
        for cell in cells:
            data = Cell.cell_data[cell]
            planes = data['z_planes']
            for plane in planes:
                xvals = np.reshape(mask_vertices[cell][plane][0], [-1, 1])
                yvals = np.reshape(mask_vertices[cell][plane][1], [-1, 1])
                zvals = np.ones([len(xvals), 1])*plane
                coords = np.concatenate((zvals, xvals, yvals), axis = 1)
                all_cell_pixels = np.concatenate((all_cell_pixels, coords), axis = 0)
        
        all_cell_pixels = all_cell_pixels[1:, :].astype(int)
            
        print('Data loaded')
        
except IOError:
    print('No saved data found, calculating mask pixels')
    c_rois = viewer.layers['Cell masks'].data
    px = pixels_in_roi.pixels_in_roi(h, w, n, c_rois)
    all_cell_pixels = px['all_pixels']
    cell_pixels_roi = px['pixels_roi']
    with open('{0}{1}{2}'.format(save_loc, sep, save_file), 'wb') as f:
        pkl.dump({'all_cell_pixels': all_cell_pixels, 'cell_pixels_roi': cell_pixels_roi}, f)
#    px = pixels_in_roi.pixels_in_roi(h, w, l_rois)
#    all_cell_pixels = px['all_pixels']


# %% Get pixels in lipofuscin ROIs

save_loc = 'G:\\Shared drives\\as_share\\HCR\\HCR_10.17'
save_file = 'HCR_10.17_S2_lipofuscin_pixels_in_cells.pkl'
    
try:
    with open('{0}{1}{2}'.format(save_loc, sep, save_file), 'rb') as f:
        mask_vertices = pkl.load(f)
        all_lipo_pixels = mask_vertices['all_lipo_pixels']
        lipo_pixels_roi = mask_vertices['lipo_pixels_roi']
        print('Data loaded')
        
except IOError:
    print('No saved data found, calculating mask pixels')
    px = pixels_in_roi.pixels_in_roi(h, w, n, l_rois)
    all_lipo_pixels = px['all_pixels']
    lipo_pixels_roi = px['pixels_roi']
    with open('{0}{1}{2}'.format(save_loc, sep, save_file), 'wb') as f:
        pkl.dump({'all_lipo_pixels': all_lipo_pixels, 'lipo_pixels_roi': lipo_pixels_roi}, f)
    
# %% Get background 
with open('{0}\\S2_background_ndnf_sst_vip.pkl'.format(save_loc), 'rb') as f:
    dict = pkl.load(f)
    
avg_bg_ndnf = dict['Ndnf'] 
avg_bg_sst = dict['Sst']  
avg_bg_vip = dict['Vip'] 
    


# %% Filter images with difference of gaussians to amplify lipofuscin-sized spots

sigma_small = 2
sigma_large = 5

t0 = time.time()    

img = im_array_ndnf
im_diff_ndnf = diff_gauss.diff_gauss(sigma_small, sigma_large, img, do_plot = 1)
t1 = time.time() - t0
print('{0} seconds'.format(int(t1)))

img = im_array_sst
im_diff_sst = diff_gauss.diff_gauss(sigma_small, sigma_large, img, do_plot = 0)
t1 = time.time() - t0
print('{0} seconds'.format(int(t1)))

img = im_array_vip
im_diff_vip = diff_gauss.diff_gauss(sigma_small, sigma_large, img, do_plot = 0)
t1 = time.time() - t0
print('{0} seconds'.format(int(t1)))

# %% 
viewer.add_image(data = im_diff_ndnf, name = 'Diff ndnf', colormap = 'magenta', blending = 'additive')
viewer.add_image(data = im_diff_sst, name = 'Diff Sst', colormap = 'green', blending = 'additive')
viewer.add_image(data = im_diff_vip, name = 'Diff Vip', colormap = 'yellow', blending = 'additive')

# %% Plot histogram of filtered images cell pixels and lipofuscin pixels

#img_dict = {'Ndnf': im_diff_ndnf, 'Sst': im_diff_sst, 'Vip': im_diff_vip}
img_dict = {'Ndnf': im_array_ndnf, 'Sst': im_array_sst, 'Vip': im_array_vip}

pixels_dict = {'Cells': all_cell_pixels, 'Lipofuscin': all_lipo_pixels}
#pixels_dict = {'Cells': all_cell_pixels}
colors_dict = {'Cells': 'b', 'Lipofuscin': 'k'}
#title = 'Diff gauss sigma = ({0}, {1})'.format(sigma_small, sigma_large)
title = 'Raw images'
save_loc = 'G:\\Shared drives\\as_share\\HCR\\HCR_10.17\\S2_scatter_plots'
save_file = '{0} hist lipo + cells.png'.format(title)
#save_file = '{0} hist cells.png'.format(title)

hist_channels.hist_channels(img_dict, pixels_dict, max_per_group= 100000, colors_dict = colors_dict, 
                            title = title, save = True, save_loc = save_loc, save_file = save_file)

del img_dict
del pixels_dict

# %% Binarize filtered images cell pixels and lipofuscin pixels

img_dict = {'Ndnf': im_diff_ndnf, 'Sst': im_diff_sst, 'Vip': im_diff_vip}
#img_dict = {'Ndnf': im_array_ndnf, 'Sst': im_array_sst, 'Vip': im_array_vip}

pixels_dict = {'Cells': all_cell_pixels, 'Lipofuscin': all_lipo_pixels}
#pixels_dict = {'Cells': all_cell_pixels}
colors_dict = {'Cells': 'b', 'Lipofuscin': 'k'}
title = 'Diff gauss sigma = ({0}, {1})'.format(sigma_small, sigma_large)
#title = 'Raw images'
save_loc = 'G:\\Shared drives\\as_share\\HCR\\HCR_10.17\\S2_scatter_plots'
save_file = '{0} hist lipo + cells.png'.format(title)
#save_file = '{0} hist cells.png'.format(title)

thresh_scale = 1

thresh = hist_channels.hist_channels(img_dict, pixels_dict, max_per_group= 100000, 
                                     do_bin = True, bin_group = 'Cells', thresh_scale = thresh_scale,
                                     colors_dict = colors_dict, 
                                     title = title, save = True, save_loc = save_loc, save_file = save_file)

del img_dict
del pixels_dict

im_bin_ndnf = np.zeros(im_diff_ndnf.shape)
im_bin_sst = np.zeros(im_diff_ndnf.shape)
im_bin_vip = np.zeros(im_diff_ndnf.shape)

ix_ndnf = np.where(im_diff_ndnf > thresh['Ndnf'])
ix_sst = np.where(im_diff_sst > thresh['Sst'])
ix_vip = np.where(im_diff_vip > thresh['Vip'])

im_bin_ndnf[ix_ndnf] = im_diff_ndnf[ix_ndnf]
im_bin_sst[ix_sst] = im_diff_sst[ix_sst]
im_bin_vip[ix_vip] = im_diff_vip[ix_vip]

# %% Find mahalanobis distance of pixels from lipofuscin cloud (in binarized images)

img_dict = {'Ndnf': im_bin_ndnf, 'Sst': im_bin_sst, 'Vip': im_bin_vip}
#img_dict = {'Ndnf': im_diff_ndnf, 'Sst': im_diff_sst, 'Vip': im_diff_vip}
#img_dict = {'Ndnf': im_array_ndnf, 'Sst': im_array_sst, 'Vip': im_array_vip}

pixels_dict = {'Cells': all_cell_pixels, 'Lipofuscin': all_lipo_pixels}

origin_group = 'Lipofuscin'
bin_group = 'Cells'

colors_dict = {'Cells': 'b', 'Lipofuscin': 'k'}

title = 'Diff gauss sigma = ({0}, {1}); binarized above {2} std from median'.format(sigma_small, sigma_large, thresh_scale)
#title = 'Raw images'
save_loc = 'G:\\Shared drives\\as_share\\HCR\\HCR_10.17\\S2_scatter_plots'
save_file = '{0} hist lipo + cells.png'.format(title)
#save_file = '{0} hist cells.png'.format(title)

output = m_dist.m_dist(img_dict, pixels_dict, origin_group, bin_group, colors_dict = colors_dict, thresh_scale = 0,
                            title = title, save = True, save_loc = save_loc, save_file = save_file)

m_dist_vals = output['m_dist']
thresh = output['thresh']

# %% Use mahalanobis distance to label all lipofuscin pixels in image (based on threshold)

lipo_pixels = all_cell_pixels[m_dist_vals['Cells'] < thresh, :]
cell_pixels = all_cell_pixels[m_dist_vals['Cells'] > thresh, :]

im_array_lipo = np.zeros(im_array_gad1.shape)
im_array_lipo[lipo_pixels[:, 0], lipo_pixels[:, 1], lipo_pixels[:, 2]] = np.ones(lipo_pixels.shape[0])

viewer.add_image(data = im_array_lipo, name = 'Detected lipofuscin', colormap = 'gray', blending = 'additive')



# %% Make scatter plots of images with different filters

img_dict = {'Ndnf': im_bin_ndnf, 'Sst': im_bin_sst, 'Vip': im_bin_vip}
#img_dict = {'Ndnf': im_diff_ndnf, 'Sst': im_diff_sst, 'Vip': im_diff_vip}
#img_dict = {'Ndnf': im_array_ndnf, 'Sst': im_array_sst, 'Vip': im_array_vip}

pixels_dict = {'Cells': all_cell_pixels, 'Lipofuscin': all_lipo_pixels}
#pixels_dict = {'Cells': all_cell_pixels}
colors_dict = {'Cells': 'b', 'Lipofuscin': 'k'}
title = 'Diff gauss sigma = ({0}, {1}); binarized above {2} std from median'.format(sigma_small, sigma_large, thresh_scale)
#title = 'Raw images'
save_loc = 'G:\\Shared drives\\as_share\\HCR\\HCR_10.17\\S2_scatter_plots'
save_file = '{0} scatter lipo + cells.png'.format(title)
#save_file = '{0} scatter cells.png'.format(title)

scatter_channels.scatter_channels(img_dict, pixels_dict, max_per_group= 100000, colors_dict = colors_dict, 
                                  title = title, save = True, save_loc = save_loc, save_file = save_file)

del img_dict
del pixels_dict

# %% Make 3D scatter plots of images with different filters

img_dict = {'Ndnf': im_bin_ndnf, 'Sst': im_bin_sst, 'Vip': im_bin_vip}
#img_dict = {'Ndnf': im_diff_ndnf, 'Sst': im_diff_sst, 'Vip': im_diff_vip}
#img_dict = {'Ndnf': im_array_ndnf, 'Sst': im_array_sst, 'Vip': im_array_vip}

pixels_dict = {'Cells': all_cell_pixels, 'Lipofuscin': all_lipo_pixels}
#pixels_dict = {'Cells': all_cell_pixels}
#pixels_dict = {'Lipofuscin': all_lipo_pixels}


#colors_dict = {'Cells': 'b', 'Lipofuscin': 'k'}
#colors_dict = {'Cells': m_dist_vals['Cells'], 'Lipofuscin': m_dist_vals['Lipofuscin']}
colors_dict = {'Cells': np.log(m_dist_vals['Cells']), 'Lipofuscin': 'r'}

title = 'Sub-sampled diff gauss sigma = ({0}, {1}) \n binarized above {2} std from median \n color = mahalanobis distance from lipofuscin pixels'.format(sigma_small, sigma_large, thresh_scale)
#title = 'Raw images'

save_loc = 'G:\\Shared drives\\as_share\\HCR\\HCR_10.17\\S2_scatter_plots'
save_file = '{0} 3D scatter lipo + cells.png'.format(title)
#save_file = '{0} 3D scatter cells.png'.format(title)

scatter_channels.scatter_channels(img_dict, pixels_dict, max_per_group= 100000, colors_dict = colors_dict,
                                  make_3D = True,
                                  title = title, save = True, save_loc = save_loc, save_file = save_file)

del img_dict
del pixels_dict

# %% Subtract background
x = 0
idx1 = 0
idx2 = 0
for cell in cells:
    data = Cell.cell_data[cell]
    planes = data['z_planes']
    
    for plane in planes:
        xvals = mask_vertices[cell][plane][0]
        idx2 = idx1 + len(xvals)
        cells_ndnf[idx1:idx2] = cells_ndnf[idx1:idx2] - avg_bg_ndnf[x]
        cells_sst[idx1:idx2] = cells_sst[idx1:idx2] - avg_bg_sst[x]
        cells_vip[idx1:idx2] = cells_vip[idx1:idx2] - avg_bg_vip[x]
        idx1 = idx2        

    x += 1
        


