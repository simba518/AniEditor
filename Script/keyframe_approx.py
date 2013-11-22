#! /usr/bin/env python

import os

homedir = os.path.expanduser('~')
project_dir = homedir+"/Workspace/AnimationEditor/";
approx_app = project_dir + "/Bin/Release/KeyfApprox";

data_dir = "Data/dinosaur/"
rest_tet = data_dir+"volmesh.abq"
rest_obj = data_dir+"render.obj"
key_obj  = data_dir+"tempt/key.obj"
interp_weights  = data_dir+"interp-weights.txt"

save_to = data_dir+"tempt/key_u"
os.system( approx_app+" "+rest_tet+" "+rest_obj+" "+key_obj+" "+interp_weights+" "+save_to+" 1e-6" )
