#! /usr/bin/env python

# this script is used to generate a series of nonlinear modes

import os

homedir = os.path.expanduser('~')
project_dir = homedir+"/Workspace/AnimationEditor/";
vol2obj = project_dir + "/Bin/Release/vol2obj";

sim_data_fold = project_dir+"/Data/tire/model/";
tetfile = sim_data_fold+"tire.abq"
objfile = sim_data_fold+"tire.obj"
saveto = sim_data_fold+"interp_weights.txt"
os.system( vol2obj + " " + tetfile + " " + objfile + " " + saveto)
