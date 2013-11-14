#! /usr/bin/env python

# this script is used to generate a series of nonlinear modes

import os

homedir = os.path.expanduser('~')
project_dir = homedir+"/Workspace/AnimationEditor/";
gennlmodes_app = project_dir + "/Script/apps/GenNLBasis";

sim_data_fold = project_dir+"/Data/mushroom/";
nlmode_gen_file = sim_data_fold + "nlmode_gen.ini";
os.system( gennlmodes_app + " " + nlmode_gen_file)
