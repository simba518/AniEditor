#! /usr/bin/env python

# Crop and resize the pictures of for the figures of the paper, then compile the
# document.

import os
from resize_images import *
from utility import *

figure_dir = "/home/simba/Workspace/"
# show an image
# image_file_name = figure_dir+"/tempt/a.jpg"
# show_image(image_file_name)

# resize images
src_dir = figure_dir+"/tempt"
target_dir = figure_dir+"/tempt/aa"
box = [0, 0, 1000, 1000]
# resize_images(src_dir, target_dir, box, "png")
resize_images(src_dir, target_dir, box, "jpg")

# compile the document
# make_doc_one_pass(latex_doc_dir)
