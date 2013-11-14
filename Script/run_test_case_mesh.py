#! /usr/bin/env python
import os
from utility import *
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------data-------------------------------------------------
init_file = "mtlopt_cen_keyW.ini"
method = "Opt_Z_AtA"

data_dir = "./Data/beam/"
app = "./Bin/Release/tempt --run_test=MtlOptTestSave "

#-----------------------------main-------------------------------------------------
os.system("cd /home/simba/Workspace/AnimationEditor/")
init_f = data_dir+init_file+" "
opt_z = data_dir+"/tempt/"+init_file+method+"OptZ.b "
opt_k = data_dir+"/tempt/"+init_file+method+"OptK.b "
cmd = app+" "+init_f+opt_z+opt_k+method
print "running: "+ cmd
os.system(cmd)
