#! /usr/bin/env python

import os
import sys

parameter_is_correct = True

rootdir = "./Bin/"
if len(sys.argv)==2:
    if sys.argv[1] == "-d":
        rootdir = rootdir + "Debug"
        print "Run Debug Test Cases\n"
    elif sys.argv[1] == "-r":
        rootdir = rootdir + "Release"
        print "Run Release Test Cases"
    else :
        print "error: the parameters for this python script is uncorrect"
        print "usage: runalltest [-r|-d]\n"
        parameter_is_correct = False

if parameter_is_correct:
    fileList = []
    for root, subFolders, files in os.walk(rootdir):
        for onefile in files:
            if (onefile[-4:]).lower() == 'test' :
                fileList.append(os.path.join(root,onefile))
    
    count = 1
    for files in fileList:
        print "\n\n----------------------------"+str(count)+"/"+str(len(fileList))+"--------------------------------"
        print "TestCase Name: " +  files
        os.system(files)
        count = count + 1
