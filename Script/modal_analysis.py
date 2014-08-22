#! /usr/bin/env python

# calculate the eigen values and eigen vectors which is used for the warping
# process.

import os

homedir = os.path.expanduser('~')
root_path = homedir + "/Workspace/AnimationEditor/"
sim_dir = root_path+"Data/flower/model/"
space = " "
vol_filename = space + sim_dir + "/flower.hinp"
fixed_nodes_file = space + sim_dir + "/con_nodes.bou"

linear_mode_num = 100
eigen_values = space+sim_dir+"/eigenvalues"+str(linear_mode_num)+".b"
eigen_vectors = space+sim_dir+"/eigenvectors"+str(linear_mode_num)+".b"

MA_app = root_path + "/Script/apps/ModalAnalysis"
os.chdir(sim_dir)
os.system(MA_app + vol_filename + eigen_values + eigen_vectors + space + str(linear_mode_num) + fixed_nodes_file)
