#! /usr/bin/env python

import os

homedir = os.path.expanduser('~')
project_dir = homedir+"/Workspace/AnimationEditor/";
approx_app = project_dir + "/Bin/Release/KeyfApprox";

data_dir = "Data/beam/"
rest_tet = data_dir+"mesh/beam.abq"
rest_obj = data_dir+"tempt/key_input.obj" 
key_obj  = data_dir+"tempt/key_output.obj"
interp_weights  = data_dir+"mesh/interp-weights.txt"

save_to = data_dir+"tempt/key_u"
os.system( approx_app+" "+rest_tet+" "+rest_obj+" "+key_obj+" "+interp_weights+" "+save_to+" 1e-3" )
