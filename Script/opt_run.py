#! /usr/bin/env python
import os
from utility import *
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------data-------------------------------------------------
init_files = [ 
               # "./Data/beam_center/mtlopt_fix.ini",
               # "./Data/beam_center/mtlopt_fix_sym.ini",
               # "./Data/beam_end/mtlopt.ini",
               # "./Data/bird/mtlopt.ini",
               "./Data/flower_one/mtlopt.ini",
               # "./Data/CharacterI/mtlopt.ini",
               ]

app = "./Bin/Release/AniEditUI "

#-----------------------------main-------------------------------------------------

os.system("cd /home/simba/Workspace/AnimationEditor/")

for inif in init_files:

    if not inif.endswith(".ini"): continue

    print "--------------------------------------------------------------------"

    model_file = str.split(inif,"/")

    cmd = app + inif + " > " + "./tempt_opt/"+model_file[-2]+"-"+model_file[-1]+".mtllog"

    print "running: "+ cmd

    os.system(cmd)
