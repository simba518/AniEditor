#! /usr/bin/env python
import os
from prepare_sim_data import *

def main():
    homedir = os.path.expanduser('~')
    script_path = homedir+"/Workspace/AnimationEditor/Script/prepare"
    data_root   = homedir+"/Workspace/AnimationEditor/Data"
    temp_fold   = data_root + "/pattern"
    obj_file    = data_root + "/flower_box_interactive.obj"

    m_SimModelGen = SimModelGen( obj_file, data_root, temp_fold, script_path )
    m_SimModelGen.mtlopt_inifile = "./Data/flower_box_interactive/mtlopt_interactive.ini"
    m_SimModelGen.startMtlOPt()

if __name__ == "__main__":
    main()
