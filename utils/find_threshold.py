"""
Created on Wed 1 Apr 2020

@author: Amrita Singh
"""

import numpy as np

def find_threshold(points, scale, type = 'Std from median'):

    center = np.median(points)
    neg_points = points[points < center]
    pos_points = neg_points + center
    points = np.concatenate([neg_points, pos_points])

    thresh = center + np.std(points)*scale

    return thresh
