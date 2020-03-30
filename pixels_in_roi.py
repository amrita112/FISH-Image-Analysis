# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 16:49:48 2020
This function finds all the 2D pixel coordinates within rois for a list of 2D rois.
Input:
    h: Image height in pixels
    w: Image width in pixels
    n: Number of planes
    rois: list of rois [vertex coordinates, 2D or 3D. If 2D but n > 1, provide plane for each roi. If 3D, shape should be (# of vertices, 3)]
    planes (if n > 1): plane for each roi
Output:
    pixels_roi: Dictionary containing list of pixel coordinates within each roi
    all_pixels: 1X2 numpy array of pixel coordinates [[x_coordinates][y_coordinates]] for all rois concatenated
@author: JRCLUST
"""
import numpy as np
import matplotlib.path as mpltpath
import time

def pixels_in_roi(h, w, n, rois, planes = None):
    
    t0 = time.time()
    n_rois = len(rois)
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
                
    t1 = time.time() - t0
    print('Calculated grid coordinates, {0} seconds'.format(np.round(t1)))
    
    if n == 1:
        
    
        pixels_roi = {}
        all_pixels = np.zeros([1, 2])
        for roi in range(n_rois):
            print('ROI {0}'.format(roi))
            
            vertices = rois[roi][:, 1:]
            path = mpltpath.Path(vertices)
            mask = path.contains_points(points)
            mask = np.reshape(mask, [h, w])
            pixels_roi[roi] = np.where(mask)
            
            n_px = len(pixels_roi[roi][0])
            xvals = np.reshape(pixels_roi[roi][0], (n_px, 1)) 
            yvals = np.reshape(pixels_roi[roi][1], (n_px, 1))
            coords = np.concatenate((xvals, yvals), axis = 1)
            all_pixels = np.concatenate((all_pixels, coords))
    
        all_pixels = all_pixels.astype(int)
        all_pixels = all_pixels[1:, :]
        
    else:
        if n > 1:
        
            pixels_roi = {}
            all_pixels = np.zeros([1, 3])
            if rois[0].shape[1] == 3:
                for roi in range(n_rois):
                    print('ROI {0}'.format(roi))
                    
                    vertices = rois[roi][:, 1:]
                    plane = rois[roi][0, 0]
                    path = mpltpath.Path(vertices)
                    mask = path.contains_points(points)
                    mask = np.reshape(mask, [h, w])
                    pixels_roi[roi] = np.where(mask)
                    
                    n_px = len(pixels_roi[roi][0])
                    xvals = np.reshape(pixels_roi[roi][0], (n_px, 1)) 
                    yvals = np.reshape(pixels_roi[roi][1], (n_px, 1))
                    coords = np.concatenate((np.ones([n_px, 1])*plane, xvals, yvals), axis = 1)
                    
                    all_pixels = np.concatenate((all_pixels, coords))
            else:
                if planes == None:
                    print('Error: please provide 3D rois or planes for rois')
                    return
                else:
                    
                    for roi in range(n_rois):
                        print('ROI {0}'.format(roi))
                        
                        vertices = rois[roi][:, 1:]
                        path = mpltpath.Path(vertices)
                        mask = path.contains_points(points)
                        mask = np.reshape(mask, [h, w])
                        pixels_roi[roi] = np.where(mask)
                        
                        n_px = len(pixels_roi[roi][0])
                        xvals = np.reshape(pixels_roi[roi][0], (n_px, 1)) 
                        yvals = np.reshape(pixels_roi[roi][1], (n_px, 1))
                        coords = np.concatenate((np.ones([1, n_px])*planes[roi], xvals, yvals), axis = 1)
                        
                        all_pixels = np.concatenate((all_pixels, coords))
        
            all_pixels = all_pixels.astype(int)
            all_pixels = all_pixels[1:, :]
        else:
            print('Error: n should be >= 1')
            return
        
    return{'pixels_roi': pixels_roi, 'all_pixels': all_pixels}