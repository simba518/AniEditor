#! /usr/bin/env python
import os
from utility import *
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------data-------------------------------------------------
init_files = ["mtlopt_cen_keyW.ini"]
test_cases = ["Opt_Z","Opt_Z_Lambda","Opt_Z_AtA"]

data_dir = "./Data/beam/"
app = "./Bin/Release/tempt --run_test=MtlOptTest/"

#-----------------------------main-------------------------------------------------
os.system("cd /home/simba/Workspace/AnimationEditor/")
for inif in init_files:
    for tc in test_cases:
        print "--------------------------------------------------------------------"
        log_f = "./tempt/"+inif+tc+".mtllog"
        cmd = app+tc+" "+data_dir+inif+" > "+log_f
        print "running: "+ cmd
        os.system(cmd)
