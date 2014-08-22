#! /usr/bin/env python

import os

homedir = os.path.expanduser('~')
project_dir = homedir+"/Workspace/AnimationEditor/";
cubature_app = project_dir + "/Bin/Release/Cubature";

sim_data_fold = project_dir+"/Data/tire/";
cubature_file = sim_data_fold + "cubature.ini";

os.system( cubature_app + " " + cubature_file)
