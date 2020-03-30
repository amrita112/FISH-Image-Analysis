# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 11:18:45 2019

@author: singha
"""

from IPython.display import display
from PIL import Image, ImageDraw
import numpy as np
import pickle as pkl
import matplotlib.path as mpltpath
import matplotlib.pyplot as plt
import cv2

## Specify images to load

folder = 'Z:\HCR\HCR_9.13\S1_546_647_dapi_514_z-stack_tiles_2019_09_16__22_04_08-Rotate 2D-01'

base_filename = 'S1_546_647_dapi_514_z-stack_tiles_2019_09_16__22_04_08-Rotate 2D-01_z'

plane_nos = range(1, 19)
n = len(plane_nos)
print('Number of planes: {0}'.format(n))

img = Image.open('{0}\{1}{2}c1_ORG.tif'.format(folder, base_filename, str(plane_nos[0]).zfill(2)))
h = img.height
w = img.width
print('Size of image in pixels: {0} X {1} X {2}'.format(n, h, w))

img_file = 'S1_546_c1_ndnf_binary'
img = Image.open('{0}\{1}.tif'.format(folder, img_file))

im_blob = cv2.imread('{0}\{1}.tif'.format(folder, img_file), 0)
# Set up the detector with default parameters.
detector = cv2.SimpleBlobDetector()
params = cv2.SimpleBlobDetector_Params()