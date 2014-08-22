#! /usr/bin/env python

import os
from prepare_sim_data import *

def main():
    homedir = os.path.expanduser('~')
    script_path = homedir+"/Workspace/AnimationEditor/Script/prepare"
    data_root   = homedir+"/Workspace/AnimationEditor/Data"
    temp_fold   = data_root + "/pattern"
    # obj_file    = data_root + "/dino_setmtl.obj"
    obj_file    = data_root + "/dino_scaled.obj"

    m_SimModelGen = SimModelGen ( obj_file, data_root, temp_fold, script_path )
    # m_SimModelGen.printAll()
    # m_SimModelGen.genSimDir()
    # m_SimModelGen.genInterpWeights()
    # m_SimModelGen.setFixedNodes()
    # m_SimModelGen.startFullStVKSimulation()
    # m_SimModelGen.startAniEditing()
    # m_SimModelGen.genNLBasis()
    # m_SimModelGen.genCubture()
    m_SimModelGen.startMtlOPt()

if __name__ == "__main__":
    main()
