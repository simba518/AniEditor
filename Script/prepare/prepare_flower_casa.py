#! /usr/bin/env python
import os
from prepare_sim_data import *

def main():
    homedir = os.path.expanduser('~')
    script_path = homedir+"/Workspace/AnimationEditor/Script/prepare"
    data_root   = homedir+"/Workspace/AnimationEditor/Data"
    temp_fold   = data_root + "/pattern"
    obj_file    = data_root + "/flower_box_casa.obj"

    m_SimModelGen = SimModelGen( obj_file, data_root, temp_fold, script_path )
    # m_SimModelGen.genNLBasis()
    # m_SimModelGen.genCubture()

    m_SimModelGen.mtlopt_inifile = "./Data/flower_box_casa/mtlopt.ini"
    m_SimModelGen.startMtlOPt()
    # m_SimModelGen.startAniEditing()

    # m_SimModelGen.aniedit_app="/home/simba/Workspace/AnimationEditing/Bin/Release/AniEditUI"
    # m_SimModelGen.aniedit_inifile = "./Data/flower_box_casa/aniedit_casa.ini"
    # m_SimModelGen.genNLBasis()
    # m_SimModelGen.startAniEditing()

if __name__ == "__main__":
    main()
