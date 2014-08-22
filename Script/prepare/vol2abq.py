#! /usr/bin/env python

import os
import string

def vol2abq(volfile_name, abqfile_name):
    """
    Convert a vol format file to abq file. The .vol format is the default output
    of the netgent when calling through commandline, and .abq is the format that
    used by my application.
    """
    # read in vol mesh
    volfile = open(volfile_name,"r")
    # read in element
    elem_list = []
    line = "begin"
    while line:
        line = volfile.readline()
        if line == "volumeelements\n":
            total_ele = string.atoi(volfile.readline())
            for i in range(1,total_ele+1):
                line = volfile.readline()
                line = line.split()
                el=str(i)+", "+line[2]+", "+line[3]+", "+line[4]+", "+line[5]+"\n"
                elem_list.append(el)

    # read in nodes
    volfile.seek(0, 0)
    node_list = []
    line = "begin"
    while line:
        line = volfile.readline()
        if line == "points\n":
            total_node = string.atoi(volfile.readline())
            print total_node
            for i in range(1,total_node+1):
                line = volfile.readline()
                line = line.split()
                node = str(i)+", "+line[0]+", "+line[1]+", "+line[2]+"\n"
                node_list.append(node)
    
    # write abq file
    abqfile = open(abqfile_name,"w")
    # write nodes
    abqfile.write("*NODE\n")
    for node in node_list:
        abqfile.write(node)
    
    # write elements
    abqfile.write("*ELEMENT, type=C3D4, ELSET=PART1\n")
    for elem in elem_list:
        abqfile.write(elem)
    abqfile.write("*ELSET,ELSET=EALL,GENERATE\n")
    abqfile.write("1,"+str(len(elem_list))+"\n")

    volfile.close()
    abqfile.close()

if __name__ == "__main__":
    volfile_name = "./Data/bird/bird247.vol"
    abqfile_name = "./Data/bird/bird247.abq"
    vol2abq(volfile_name, abqfile_name)
