#! /usr/bin/env python

# this script is used to generate a series of nonlinear modes

import os

homedir = os.path.expanduser('~')
project_dir = " "+homedir+"/Workspace/AnimationEditor/";
convert_app = project_dir + "/Bin/Release/key2par ";

sim_data_fold = project_dir+"/Data/beam/W80_B40_C97_FixEnd/"

tetmesh_file = project_dir+"/Data/beam/mesh/beam.abq "
keyframes_file = sim_data_fold+"keyU4.b " 
cub_points_file = sim_data_fold+"cubaturePoints.txt " 
cub_weights_file = sim_data_fold+"cubatureWeights.b " 
num_desired_points = " 100 "
keyframes = " 0 10 20 30 "

os.system( convert_app+tetmesh_file+keyframes_file+cub_points_file+cub_weights_file+keyframes+num_desired_points )
