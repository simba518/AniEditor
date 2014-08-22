#! /usr/bin/env python
import os
from utility import *
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------data-------------------------------------------------
init_file = "mtlopt.ini"
method = "Opt_Z_Lambda"

data_dir = "./Data/beam_center/"
app = "./Bin/Release/MtlOptTest --run_test=MtlOptTestSave "

#-----------------------------main-------------------------------------------------
os.system("cd /home/simba/Workspace/AnimationEditor/")
init_f = data_dir+init_file+" "
opt_z = data_dir+"/tempt/"+init_file+method+"OptZ.b "
opt_k = data_dir+"/tempt/"+init_file+method+"OptK.b "
opt_s = data_dir+"/tempt/"+init_file+method+"OptS.b "
cmd = app+" "+init_f+opt_z+opt_k+opt_s+method
print "running: "+ cmd
os.system(cmd)
