#! /usr/bin/env python

# this script is used to generate a series of nonlinear modes

import os

homedir = os.path.expanduser('~')
project_dir = homedir+"/Workspace/AnimationEditor/";
cubature_app = project_dir + "/Script/apps/Cubature";

sim_data_fold = project_dir+"/Data/beam/";

# warpbasis_gen_file = sim_data_fold + "warpbasis.ini";
# nlmode_gen_file = sim_data_fold + "nlmode-gen.ini";
cubature_file = sim_data_fold + "cubature.ini";

os.system( cubature_app + " " + cubature_file)
