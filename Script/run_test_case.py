#! /usr/bin/env python
import os
from utility import *
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------data-------------------------------------------------
init_files = ["mtlopt2.ini"]

test_cases = ["Opt_Z_AtA"]
data_dir = "./Data/mushroom/"
app = "./Bin/Release/MtlOptTest --run_test=MtlOptTest/"

#-----------------------------main-------------------------------------------------
os.system("cd /home/simba/Workspace/AnimationEditor/")
for inif in init_files:
    if not inif.endswith(".ini"): continue
    for tc in test_cases:
        print "--------------------------------------------------------------------"
        log_f = "./tempt/"+inif+tc+".mtllog"
        cmd = app+tc+" "+data_dir+inif+" > "+log_f
        print "running: "+ cmd
        os.system(cmd)
