#! /usr/bin/env python
import os
from utility import *
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------data-------------------------------------------------
report_tex = "./tempt/report.tex"

#-----------------------------functions--------------------------------------------
def test_error_check(log_file):
    succ = True
    f = open(log_file)
    for line in f:
        if line.find("failure") >= 0:
            succ = False
            break
    return succ

def model_info(log_file):
    print log_file
    data = []
    data += grepNumberWithKey(log_file,"number of nodes:")
    data += grepNumberWithKey(log_file,"number of tetrahedrons:")
    data += grepNumberWithKey(log_file,"number of selected modes:")
    data.append(len(grepNumberWithKey(log_file,"OUTTER ITERATION")))
    data += grepNumberWithKey(log_file,"time step:")
    data += grepNumberWithKey(log_file,"alphak:")
    data += grepNumberWithKey(log_file,"alpham:")
    data += grepStrWithKey(log_file,"init file:")
    data += grepStrWithKey(log_file,"test opt method:")
    data.append(map(int,grepNumberWithKey(log_file,"keyframe id:")))
    data.append(str('%.5g'%(grepNumberWithKey(log_file,"Objective...............:")[-1])))
    data += grepNumberWithKey(log_file,"number of modes:")
    data += grepNumberWithKey(log_file,"partial penalty:")
    data.append(str(('%.5g'%grepNumberWithKey(log_file,"S.norm() =")[-1])))
    return data

def mtl_info(log_file):
    lambda_diag_d = ""
    ek0 = grepNumberWithKey(log_file,"initial values:")
    lambda_diag_d += str(map(lambda n: float('%.5g'%n), ek0))[1:-2]+"\\\ \hline\n&"
    r = len(ek0)

    ek_all = grepNumberWithKey(log_file,"eigen(K):")
    ek = []
    if(len(ek_all)>=r):
        ek += ek_all[0:r]
    if(len(ek_all)>35*r):
        ek += ek_all[-35*r:len(ek_all)]
    
    for i in range(0,len(ek)):
        lambda_diag_d += str('%.5g'%ek[i])
        if (i+1)%len(ek0) == 0:
            lambda_diag_d += "\\\ \hline\n"
            if i != len(ek)-1:
                lambda_diag_d += "&"
        else:
            lambda_diag_d += ", "
    return lambda_diag_d

def save_fig(y,fname,f_label,style='r'):
    if len(y) > 0:
        plt.semilogy()
        plt.xlim(-len(y)/20,len(y)+len(y)/20)
        plt.grid(True)
        plt.plot(range(0,len(y)),y,style,label=f_label)
    plt.savefig(fname)
    plt.clf()

def save_color_segements(y,end_of_segments,fname,colors):
    if len(y) > 0:
        plt.semilogy()
        plt.xlim(-len(y)/20,len(y)+len(y)/20)
        plt.grid(True)
        for i in range(1,len(end_of_segments)-1):
            p0 = end_of_segments[i-1]
            p1 = end_of_segments[i]
            plt.plot(range(p0,p1),y[p0:p1],colors[(i-1)%len(colors)])

    plt.savefig(fname)
    plt.clf()
    
def opt_fig(log_file):
    lof = os.path.dirname(log_file)+"/"+os.path.basename(log_file).replace(".","-")
    fig_file_name = [lof+'-inner_py.png',lof+'-outer_py.png']
    inner_y = []
    outer_y = []
    end_of_segments = [0]
    f = open(log_file)
    for line in f:
        if line.find("iter    objective") >= 0:
            for line in f:
                if len(line) <= 4: break
                if line.find("iter") < 0:
                    inner_y.append( float(line.split()[1]) )
        elif line.find("Objective...............:") >= 0:
            outer_y.append( float(line.split()[2]) )
            end_of_segments.append(len(inner_y)-1)

    # save_fig(inner_y,fig_file_name[0],'inner iteration','r')
    colors = ['r','g','b']
    save_color_segements(inner_y,end_of_segments,fig_file_name[0],colors)
    save_fig(outer_y,fig_file_name[1],'outer iteration','-b')
    return fig_file_name

def z_fig(log_file):
    all_curve_new_z = grepNumberWithKey(log_file,"curve_new_z:")
    keyframes = map(int,grepNumberWithKey(log_file,"keyframe id:"))
    if len(all_curve_new_z) > 1:
        T = int(all_curve_new_z[0])
        r = int(all_curve_new_z[1])
        curve_new_z = all_curve_new_z[len(all_curve_new_z)-T*r:len(all_curve_new_z)]
        for i in range(0,r):
            zi = curve_new_z[i*T:(i+1)*T]
            plt.grid(True)
            plt.plot(range(0,len(zi)),zi,label='m'+str(i))
            if len(keyframes) > 0:
                kzf = []
                for j in keyframes:
                    kzf.append(int(zi[j]))
                plt.plot(keyframes,kzf,'ro')

    zf=os.path.dirname(log_file)+"/"+os.path.basename(log_file).replace(".","-")+'-z.png'
    plt.legend(bbox_to_anchor=(0,1.01,1,1), loc=3, ncol=5, mode="expand",borderaxespad=0.)
    plt.savefig(zf)
    plt.clf()
    return zf

#-----------------------------main-------------------------------------------------
os.system("cd /home/simba/Workspace/AnimationEditor/")
tex_str = open("./Script/patterns/test_report_head.tex").read()
log_fs = os.listdir("./tempt")

for log_f in log_fs:
    if not log_f.endswith(".mtllog"): continue
    # if log_f.find("Z_Lambda_u") < 0: continue
    # if log_f.find("mtlopt.ini") < 0: continue
    log_f = "./tempt/"+log_f
    mode_data = model_info(log_f)
    mtl_data = mtl_info(log_f)
    energy_f = opt_fig(log_f)
    z_f = z_fig(log_f)
    tempt = report_tex+"tmpt"
    os.system("cp "+"./Script/patterns/test_report_data.tex "+tempt)
    changeElements(tempt,"#number_of_nodes#",str(int(mode_data[0])))
    changeElements(tempt,"#number_of_tets#",str(int(mode_data[1])))
    changeElements(tempt,"#number_of_outer_its#",str(int(mode_data[3])))
    changeElements(tempt,"#time_step#",str(str(mode_data[4])))
    changeElements(tempt,"#alphak#",str(str(mode_data[5])))
    changeElements(tempt,"#alpham#",str(str(mode_data[6])))
    changeElements(tempt,"#init_file_name#",str(str(mode_data[7])).replace("_","\\_"))
    changeElements(tempt,"#method#",str(str(mode_data[8])).replace("_","\\_"))
    changeElements(tempt,"#keyframes#",str(mode_data[9]))
    changeElements(tempt,"#optimal fun#",str(mode_data[10]))
    changeElements(tempt,"#number_of_modes#",str(int(mode_data[11]))+","+str(int(mode_data[2])))
    changeElements(tempt,"#partial_penalty#",str(mode_data[12]))
    if(len(mode_data)>=14):
        changeElements(tempt,"#S.norm()#",str(mode_data[13]))
    else:
        changeElements(tempt,"#S.norm()#","")

    changeElements(tempt,"#lambda_diag_d#",mtl_data)
    changeElements(tempt,"#inner_it#",energy_f[0].replace("_","\\string_"))
    changeElements(tempt,"#outer_it#",energy_f[1].replace("_","\\string_"))
    changeElements(tempt,"#curve_z#",z_f.replace("_","\\string_"))
    tex_str += open(tempt).read()

open(report_tex,'w').write(tex_str+"\n\end{document}")
os.system("pdflatex "+" -output-directory="+os.path.dirname(report_tex)+" "+report_tex)
