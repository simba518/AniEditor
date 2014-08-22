#! /usr/bin/env python

import os
from vol2abq import *
from binarySTL2asciiSTL import *

# change elements of a file
def changeElements(filename,ele_name,ele_value):
    fread = open(filename,"r")
    contents = fread.read()
    contents = contents.replace(ele_name,ele_value)
    fread.close()
    
    fwrite = open(filename,"w")
    fwrite.write(contents)
    fwrite.close()

class SimModelGen:
    
    def __init__( self, _objfile_path, _sim_data_fold_root, _pattern_fold, _script_path ):
        self.initFileAndDir( _objfile_path, _sim_data_fold_root, _pattern_fold, _script_path )

    def initFileAndDir( self, _objfile, _sim_data_fold_root, _pattern_fold, _script_path ):

        self.objfile          = _objfile
        self.objfile_path     = os.path.dirname(self.objfile)
        self.objfile_name_ext = os.path.basename(self.objfile)
        self.objfile_name     = os.path.splitext(self.objfile_name_ext)[0]
        self.backup_path      = self.objfile_path + "/backup"

        self.sim_data_fold_root = _sim_data_fold_root
        self.sim_data_fold      = self.sim_data_fold_root + "/" + self.objfile_name
        self.rendermesh_fold    = self.sim_data_fold   + "/model"
        self.rendermesh         = self.rendermesh_fold + "/mesh.obj"
        self.simp_mesh          = self.rendermesh_fold + "/mesh_simplified.obj"
        self.inflated_mesh      = self.rendermesh_fold + "/mesh_inflated.obj"
        self.bin_stl_mesh       = self.rendermesh_fold + "/mesh.stl"
        self.ascii_stl_mesh     = self.bin_stl_mesh
        self.vol_mesh           = self.sim_data_fold   + "/model/mesh.vol"
        self.abq_mesh           = self.sim_data_fold   + "/model/mesh.abq"
        self.interp_weights     = self.sim_data_fold   + "/model/interp_weights.txt"

        self.pattern_fold    = _pattern_fold
        self.script_path     = _script_path
        self.unitize_sp      = self.script_path   + "/unitize.mlx"
        self.simplify_sp     = self.script_path   + "/simplify.mlx"
        self.inflate_app     = self.script_path   + "/inflator"
        self.simulation_app  = self.script_path   + "/set_elastic_material"
        self.gennlmodes_app  = self.script_path   + "/GenNLBasis"
        self.sim_inifile     = self.sim_data_fold + "/simu_ma.ini"
        self.nlmode_gen_file = self.sim_data_fold + "/nlmodegen.ini"
        self.cubature_file   = self.sim_data_fold + "/cubature.ini"
        self.aniedit_inifile = self.sim_data_fold + "/aniedit.ini"
        self.mtlopt_inifile  = self.sim_data_fold + "/mtlopt.ini"
        self.full_stvk_sim_inifile = self.sim_data_fold + "/simu_stvk_full.ini"

    def genSimDir(self):
        sim_data_fold = os.path.expanduser(self.sim_data_fold)
        if os.path.exists(sim_data_fold):
            # os.system("rm "+sim_data_fold+" -f -r")
            print "error: the target fold is existed: "+sim_data_fold
            return

        os.system("cp "+self.pattern_fold+" "+self.sim_data_fold+" -r")

        backup_path = os.path.expanduser(self.backup_path)
        if not os.path.exists(backup_path):
            os.system("mkdir " + backup_path)
            os.system("cp " + self.objfile_path+"/*.* "+backup_path + " -f")
        
        os.system( "cp "+ self.objfile_path + "/*.* " + self.rendermesh_fold )
        os.system( "rm "+ self.rendermesh_fold + "/*.obj" )
        os.system( "cp "+ self.objfile + " " + self.rendermesh )

        init_file_dir = self.sim_data_fold + "/"
        for root, sub_folders, files in os.walk(sim_data_fold):
            for one_file in files:
                filename = os.path.join(root,one_file)
                changeElements(filename, "#model_name#", self.objfile_name)

    def chanageModelName(self):
        for root, sub_folders, files in os.walk(self.sim_data_fold):
            for one_file in files:
                filename = os.path.join(root,one_file)
                changeElements(filename, "#model_name#", self.objfile_name)


    def unitize(self):
        meshlab_cmd = " -i "+self.rendermesh+" -o "+self.rendermesh+" -s "+self.unitize_sp
        os.system("meshlabserver " + meshlab_cmd)


    def simplify(self):
        meshlab_cmd=" -i "+self.rendermesh+" -o "+self.simp_mesh+" -s "+self.simplify_sp
        os.system("meshlabserver " + meshlab_cmd)


    def inflate(self):
        if not os.path.exists(self.simp_mesh):
            self.simp_mesh = self.rendermesh
        os.system( self.inflate_app+" "+ self.simp_mesh +" "+self.inflated_mesh+" "+" 3 100 ")


    def remesh(self):
        if not os.path.exists(self.inflated_mesh):
            self.inflated_mesh = self.rendermesh
        os.system( self.remesh_app+" "+ self.inflated_mesh)

    
    def genStl(self):
        if not os.path.exists(self.inflated_mesh):
            self.inflated_mesh = self.rendermesh
        meshlab_cmd=" -i "+self.inflated_mesh+" -o "+self.bin_stl_mesh+" -s "+self.simplify_sp
        os.system("meshlabserver " + meshlab_cmd)
        binarySTL2asciiSTL( self.bin_stl_mesh, self.ascii_stl_mesh )


    def genVolMesh(self, showGui=True):
        net_cmd = "netgen -geofile=" + self.ascii_stl_mesh + " -meshfile=" + self.vol_mesh
        if not showGui:
            net_cmd = net_cmd + " -V -batchmode -coarse"
        os.system("export NETGENDIR=/usr/share/netgen/")
        os.system(net_cmd)

    def genAbqMesh(self):
        vol2abq( self.vol_mesh , self.abq_mesh)

    def genInterpWeights(self):
        app = self.vol2obj_app
        abq = self.abq_mesh
        obj = self.rendermesh
        interp = self.interp_weights
        os.system(app+" "+abq+" "+obj+" "+interp)


    def setFixedNodes(self):
        os.system(self.simulation_app + " " + self.sim_inifile + " auto_save")
        # os.system(self.simulation_app + " " + self.sim_inifile)


    def genNLBasis(self):
        app = self.gennlmodes_app
        inifile = self.nlmode_gen_file
        os.system( app + " " + inifile)


    def genCubture(self):
        app = self.cubature_app
        inifile = self.cubature_file
        os.system( app + " " + inifile)


    def startFullStVKSimulation(self):
        os.system(self.full_stvk_simulation_app + " " + self.full_stvk_sim_inifile)


    def startAniEditing(self):
        os.system(self.aniedit_app + " " + self.aniedit_inifile)


    def startMtlOPt(self):
        inif = self.mtlopt_inifile
        model_file = str.split(inif,"/")
        cmd = self.mtlopt_app+ " " + inif + " > " + "./tempt_opt/"+model_file[-2]+"-"+model_file[-1]+".mtllog"
        os.system(cmd)

    def printAll(self):
        print  "objfile            = " + self.objfile
        print  "objfile_path       = " + self.objfile_path
        print  "objfile_name_ext   = " + self.objfile_name_ext
        print  "objfile_name       = " + self.objfile_name
        print  "backup_path        = " + self.backup_path

        print  "sim_data_fold_root = " + self.sim_data_fold_root
        print  "sim_data_fold      = " + self.sim_data_fold
        print  "rendermesh_fold    = " + self.rendermesh_fold
        print  "rendermesh         = " + self.rendermesh
        print  "simp_mesh          = " + self.simp_mesh
        print  "inflated_mesh      = " + self.inflated_mesh
        print  "bin_stl_mesh       = " + self.bin_stl_mesh
        print  "ascii_stl_mesh     = " + self.ascii_stl_mesh
        print  "vol_mesh           = " + self.vol_mesh
        print  "abq_mesh           = " + self.abq_mesh
        print  "interp_weights     = " + self.interp_weights

        print  "pattern_fold       = " + self.pattern_fold
        
        print  "script_path        = " + self.script_path
        print  "unitize_sp         = " + self.unitize_sp
        print  "simplify_sp        = " + self.simplify_sp
        print  "inflate_app        = " + self.inflate_app
        print  "vol2obj_app        = " + self.vol2obj_app
        print  "simulation_app     = " + self.simulation_app
        print  "sim_inifile        = " + self.sim_inifile
        print  "gennlmodes_app     = " + self.gennlmodes_app
        print  "nlmode_gen_file    = " + self.nlmode_gen_file


    # file names and directories
    objfile          = ""
    objfile_path     = ""
    objfile_name_ext = ""
    objfile_name     = ""
    backup_path      = ""

    sim_data_fold_root = ""
    sim_data_fold      = ""
    rendermesh_fold    = ""
    
    rendermesh     = ""
    simp_mesh      = ""
    inflated_mesh  = ""
    bin_stl_mesh   = ""
    ascii_stl_mesh = ""
    vol_mesh       = ""
    abq_mesh       = ""
    interp_weights = ""
    
    pattern_fold   = ""
            
    # scripts and applications
    script_path     = ""
    unitize_sp      = "unitize.mlx"
    simplify_sp     = "simplify.mlx"
    inflate_app     = "inflator"
    vol2obj_app     = "./Bin/Release/vol2obj"
    simulation_app  = "set_elastic_material"
    sim_inifile     = ""
    full_stvk_simulation_app  = "./Bin/Release/StVKSimulator"
    full_stvk_sim_inifile     = ""
    gennlmodes_app  = "GenNLBasis"
    nlmode_gen_file = ""
    cubature_app    = "./Bin/Release/Cubature";
    cubature_file   = ""
    aniedit_app = "./Bin/Release/AniEditUI"
    aniedit_inifile = ""
    mtlopt_app = "./Bin/Release/AniEditUI"
    mtlopt_inifile = ""
    remesh_app = "/home/simba/Workspace/OthersProjects/ReMESH_2.1_linux32/remesh21"
