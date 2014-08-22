#! /usr/bin/env python

import os
from prepare_sim_data import *

def main():
    homedir = os.path.expanduser('~')
    script_path = homedir+"/Workspace/AnimationEditor/Script/prepare"
    data_root   = homedir+"/Workspace/AnimationEditor/Data"
    temp_fold   = data_root + "/pattern"
    objfiles = ["bird_small","flower_box","beam_fine","dino_uniform"]
    # objfiles = ["dino_uniform"]
    
    for objf in objfiles:
        obj_file    = data_root + "/" + objf + ".obj"
        m_SimModelGen = SimModelGen ( obj_file, data_root, temp_fold, script_path )
        # m_SimModelGen.printAll()
        # m_SimModelGen.genSimDir()
        # m_SimModelGen.unitize()
        # m_SimModelGen.simplify()
        # m_SimModelGen.inflate()
        # m_SimModelGen.remesh()
        # m_SimModelGen.genStl()
        # m_SimModelGen.genVolMesh(True)
        # m_SimModelGen.genAbqMesh()
        # m_SimModelGen.genInterpWeights()
        # m_SimModelGen.setFixedNodes()
        # m_SimModelGen.startFullStVKSimulation()
        # m_SimModelGen.startAniEditing()
        # m_SimModelGen.genNLBasis()
        # m_SimModelGen.genCubture()
        m_SimModelGen.startMtlOPt()

if __name__ == "__main__":
    main()
