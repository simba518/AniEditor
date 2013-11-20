#! /usr/bin/env python
import os
from utility import *
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------data-------------------------------------------------
# init_files = ["mtlopt_cen_keyW5.ini","mtlopt_cen_keyW_no_up5.ini",
#               "mtlopt_cen_keyW3.ini","mtlopt_cen_keyW_no_up3.ini",
#               "mtlopt_cen_keyW.ini","mtlopt_cen_keyW_no_up.ini"]

# init_files = os.listdir(data_dir)
init_files = ["mtlopt_cen_keyW3_4_new.ini",
              "mtlopt_cen_keyW3_new.ini",
              # ,"mtlopt_cen_keyW3_unp.ini",
              # "mtlopt_cen_keyW6_new.ini","mtlopt_cen_keyW6_unp.ini",
              # "mtlopt_cen_keyW8_new.ini","mtlopt_cen_keyW8_unp.ini",
              # "mtlopt_cen_keyW10_new.ini","mtlopt_cen_keyW10_unp.ini"
              ]

test_cases = [# "Opt_Z","Opt_Z_Lambda",
              "Opt_Z_AtA"]
data_dir = "./Data/beam/KeyW_new/"
app = "./Bin/Release/TemptTest --run_test=MtlOptTest/"

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
