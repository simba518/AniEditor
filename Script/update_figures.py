#! /usr/bin/env python

# Crop and resize the pictures of for the figures of the paper, then compile the
# document.

import os
from resize_images import *
from data_dir import *
from utility import *

figure_dir = "~/Workspace/AnimationEditor/Doc/"
# show an image
# image_file_name = figure_dir+"/tempt/tight_pos_con_linear.jpg"
# show_image(image_file_name)

# resize images
src_dir = figure_dir+"/tempt"
target_dir = figure_dir+"/report"
box = [380, 230, 1400, 820]
resize_images(src_dir, target_dir, box, "png")
resize_images(src_dir, target_dir, box, "jpg")

# compile the document
# make_doc_one_pass(latex_doc_dir)
