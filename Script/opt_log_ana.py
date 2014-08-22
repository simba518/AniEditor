#! /usr/bin/env python
import os
from utility import *
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------functions--------------------------------------------

def grepInt(log_file,keyname,default=""):
    data = grepNumberWithKey(log_file,keyname)
    if len(data) <= 0:
        return default
    return str(int(data[0]))


def grepFloat(log_file,keyname,default=""):
    data = grepNumberWithKey(log_file,keyname)
    if len(data) <= 0:
        return default
    return str('%.5g'%(data[0]))


def grepStr(log_file,keyname,default=""):
    data = grepStrWithKey(log_file,keyname)
    if len(data) <= 0:
        return default
    return str(data[0])

def grepFloatList(log_file,keyname,default=""):
    data = grepNumberWithKey(log_file,keyname)
    if len(data) <= 0:
        return default
    d = ""
    for var in data:
        d = d +str('%.3g'%(var))+", "
    return d

def count(log_file,key):
    data_file = open(log_file)
    data_lines = data_file.readlines()
    num = 0
    for line in data_lines:
        if line.find(key) >= 0:
            num = num+1
    return num

def sum(log_file,key):
    s = 0
    data = grepNumberWithKey(log_file,key)
    for var in data:
        s = s + var
    return str(s)

def sumInt(log_file,key):
    s = 0
    data = grepNumberWithKey(log_file,key)
    for var in data:
        s = s + int(var)
    return str(s)

def z_fig(log_file):

    data_z = ""
    all_curve_new_z = grepNumberWithKey(log_file,"curve_new_z:")
    if len(all_curve_new_z) > 1:
        T = int(all_curve_new_z[0])
        r = int(all_curve_new_z[1])
        curve_new_z = all_curve_new_z[2:-1]
        for i in range(0,7):
            zi = curve_new_z[i*T:(i+1)*T]
            data_z += str(zi)
            plt.grid(True)
            plt.plot(range(0,len(zi)),zi,label='Z$_'+str(i)+'$',linewidth=3)
            plt.xlabel("Frame number",fontsize=14)
            plt.ylabel("Z",fontsize=14)

    zf=os.path.dirname(log_file)+"/"+os.path.basename(log_file).replace(".","-")+'-z.png'
    if len(grepNumberWithKey(log_file,"optimize mtl:")) > 0:
        zf = os.path.dirname(log_file)+"/"+os.path.basename(log_file).replace(".","-")+'-z_mtlopt.png'
    legend = plt.legend(bbox_to_anchor=(0,1.01,1,20), loc=8, ncol=8, mode="expand",borderaxespad=-0.5,labelspacing=-0.2)
    for label in legend.get_texts():
        label.set_fontsize(14)
    plt.savefig(zf)
    plt.clf()

    data_z = data_z.replace("[","").replace("]","\n\n")
    df = open(zf+".txt","w")
    df.write(data_z)
    return zf


def save_fig(y,fname,f_label,style='r'):
    if len(y) > 0:
        plt.semilogy()
        plt.xlim(-len(y)/20,len(y)+len(y)/20)
        plt.grid(True)
        plt.plot(range(0,len(y)),y,style,label=f_label,linewidth=5)
    plt.savefig(fname, dpi=100)
    plt.clf()



def save_color_segements(y,end_of_segments,fname,colors):
    if len(y) > 0:
        plt.semilogy()
        plt.xlim(-len(y)/20,len(y)+len(y)/20)
        for i in range(1,len(end_of_segments)-1):
            p0 = end_of_segments[i-1]
            p1 = end_of_segments[i]+1
            plt.plot(range(p0,p1),y[p0:p1],colors[i-1],linewidth=3)

    plt.savefig(fname, dpi=100)
    plt.clf()


def controlForces(log_file):
    saveto = log_file+"-Ef.png"
    if len(grepNumberWithKey(log_file,"optimize mtl:")) > 0:
        saveto = log_file+"-Ef-mtlopt.png"

    Ef = grepNumberWithKey(log_file,"Ef = :")
    plt.grid(True)
    if len(Ef) > 0:
        # plt.plot(range(0,len(Ef)),Ef,label="Ef")
        plt.plot(range(0,len(Ef)/2),Ef[0:len(Ef)/2],label="Ef",linewidth=5)
        plt.savefig(saveto)
        plt.clf()

    Effile = open(saveto+".txt",'w')
    Effile.write(str(Ef).replace(","," ").replace("[","").replace("]",""))


def keyframeDiff(log_file):
    saveto = log_file+"-Ek.png"
    if len(grepNumberWithKey(log_file,"optimize mtl:")) > 0:
        saveto = log_file+"-Ek-mtlopt.png"

    Ef = grepNumberWithKey(log_file,"Ek = :")
    plt.grid(True)
    if len(Ef) > 0:
        plt.plot(range(0,len(Ef)),Ef,label="Ek",linewidth=5)
        # plt.plot(range(0,len(Ef)/2),Ef[0:len(Ef)/2],label="Ek",linewidth=5)
        plt.savefig(saveto)
        plt.clf()

    Effile = open(saveto+".txt",'w')
    Effile.write(str(Ef).replace(","," ").replace("[","").replace("]",""))


def Ef_add_Ek(log_file):
    saveto = log_file+"-Ek_add_Ec.png"
    if len(grepNumberWithKey(log_file,"optimize mtl:")) > 0:
        saveto = log_file+"-Ek_add_Ec-mtlopt.png"

    Ef = grepNumberWithKey(log_file,"Ef = :")
    Ek = grepNumberWithKey(log_file,"Ek = :")
    Ekf = Ek
    for i in range(0,min(len(Ef),len(Ekf))):
        Ekf[i] += Ef[i]

    plt.grid(True)
    if len(Ef) > 0:
        plt.plot(range(0,len(Ekf)),Ekf,label="Ek+Ef",linewidth=5)
        # plt.plot(range(0,len(Ef)/2),Ef[0:len(Ef)/2],label="Ek",linewidth=5)
        plt.savefig(saveto)
        plt.clf()

    Effile = open(saveto+".txt",'w')
    Effile.write(str(Ef).replace(","," ").replace("[","").replace("]",""))


def opt_fig(log_file):

    lof = os.path.dirname(log_file)+"/"+os.path.basename(log_file).replace(".","-")
    fig_file_name = [lof+'-inner_py.png',lof+'-outer_py.png']

    inner_y = []
    end_of_segments = [0]
    colors = []

    f = open(log_file)

    for line in f:
        if line.find("inner-it-fun(Z)") >= 0:
            inner_y.append( float(line.split()[2]) )
            for line in f:
                if line.find("outer-it-fun(Z)") >= 0: 
                    outer_y.append(float(line.split()[3]))
                    end_of_segments.append(len(inner_y)-1)
                    colors.append('r')
                    break
                if line.find("inner-it-fun(Z)") >= 0:
                    if len(inner_y) > 0 and float(line.split()[2]) > inner_y[-1]:
                        inner_y.append( inner_y[-1] )
                    else:
                        inner_y.append( float(line.split()[2]) )

        if line.find("inner-it-fun-grad(S)") >= 0:
            inner_y.append( float(line.split()[3]) )
            for line in f:
                if line.find("outer-it-fun(S)") >= 0:
                    outer_y.append( float(line.split()[3]) )
                    end_of_segments.append(len(inner_y)-1)
                    colors.append('g')
                    break
                if line.find("inner-it-fun-grad(S)") >= 0:
                    if len(inner_y) > 0 and float(line.split()[3]) > inner_y[-1]:
                        inner_y.append( inner_y[-1] )
                    else:
                        inner_y.append( float(line.split()[3]) )

        if line.find("inner-it-fun(Lambda)") >= 0:
            inner_y.append( float(line.split()[3]) )
            for line in f:
                if line.find("outer-it-fun(Lambda)") >= 0:
                    outer_y.append( float(line.split()[3]) )
                    end_of_segments.append(len(inner_y)-1)
                    colors.append('b')
                    break
                if line.find("inner-it-fun(Lambda)") >= 0:
                    if len(inner_y) > 0 and float(line.split()[3]) > inner_y[-1]:
                        inner_y.append( inner_y[-1] )
                    else:
                        inner_y.append( float(line.split()[3]) )

    if len(inner_y) > 0:
        save_color_segements(inner_y,end_of_segments,fig_file_name[0],colors)
        outer_y[0] = inner_y[0]
        save_fig(outer_y,fig_file_name[1],'outer iteration','-b')

    inner_f = open(fig_file_name[0]+".txt",'w')
    inner_f.write(str(inner_y).replace(","," ").replace("[","").replace("]",""))
    inner_f = open(fig_file_name[0]+"_ends.txt",'w')
    inner_f.write(str(end_of_segments).replace(","," ").replace("[","").replace("]",""))

    return fig_file_name

    

def opt_fig(log_file):

    lof = os.path.dirname(log_file)+"/"+os.path.basename(log_file).replace(".","-")
    fig_file_name = [lof+'-inner_py.png',lof+'-outer_py.png']

    outer_y = [0]
    inner_y = []
    end_of_segments = [0]
    colors = []

    f = open(log_file)

    for line in f:
        if line.find("inner-it-fun(Z)") >= 0:
            inner_y.append( float(line.split()[2]) )
            for line in f:
                if line.find("outer-it-fun(Z)") >= 0: 
                    outer_y.append(float(line.split()[3]))
                    end_of_segments.append(len(inner_y)-1)
                    colors.append('r')
                    break
                if line.find("inner-it-fun(Z)") >= 0:
                    if len(inner_y) > 0 and float(line.split()[2]) > inner_y[-1]:
                        inner_y.append( inner_y[-1] )
                    else:
                        inner_y.append( float(line.split()[2]) )

        if line.find("inner-it-fun-grad(S)") >= 0:
            inner_y.append( float(line.split()[3]) )
            for line in f:
                if line.find("outer-it-fun(S)") >= 0:
                    outer_y.append( float(line.split()[3]) )
                    end_of_segments.append(len(inner_y)-1)
                    colors.append('g')
                    break
                if line.find("inner-it-fun-grad(S)") >= 0:
                    if len(inner_y) > 0 and float(line.split()[3]) > inner_y[-1]:
                        inner_y.append( inner_y[-1] )
                    else:
                        inner_y.append( float(line.split()[3]) )

        if line.find("inner-it-fun(Lambda)") >= 0:
            inner_y.append( float(line.split()[3]) )
            for line in f:
                if line.find("outer-it-fun(Lambda)") >= 0:
                    outer_y.append( float(line.split()[3]) )
                    end_of_segments.append(len(inner_y)-1)
                    colors.append('b')
                    break
                if line.find("inner-it-fun(Lambda)") >= 0:
                    if len(inner_y) > 0 and float(line.split()[3]) > inner_y[-1]:
                        inner_y.append( inner_y[-1] )
                    else:
                        inner_y.append( float(line.split()[3]) )

        # if line.find("iter    objective    inf_pr") >= 0:
        #     for line in f:
        #         if len(line) <= 10 or line.find("Number of Iterations....")>=0:
        #             end_of_segments.append(len(inner_y))
        #             colors.append('b')
        #             break
        #         if (line[0] != 'i'):
        #             inner_y.append( float(line.split()[1]) )

    if len(inner_y) > 0:
        save_color_segements(inner_y,end_of_segments,fig_file_name[0],colors)
        outer_y[0] = inner_y[0]
        save_fig(outer_y,fig_file_name[1],'outer iteration','-b')

    inner_f = open(fig_file_name[0]+".txt",'w')
    inner_f.write(str(inner_y).replace(","," ").replace("[","").replace("]",""))
    inner_f = open(fig_file_name[0]+"_ends.txt",'w')
    inner_f.write(str(end_of_segments).replace(","," ").replace("[","").replace("]",""))

    return fig_file_name

#-----------------------------main------------------------------------------------
report_tex = "./tempt_opt/report.tex"

log_fs = os.listdir("./tempt_opt")

tex_str = open("./Script/patterns/opt_report_head.tex").read()

for log_f in log_fs:

    if not log_f.endswith(".mtllog"): continue
    if not log_f.find("flower_box") >= 0: continue

    tempt = report_tex+"tmpt"
    os.system("cp "+"./Script/patterns/opt_report_data.tex "+tempt)

    log_f = "./tempt_opt/"+log_f

    z_f = z_fig(log_f)
    changeElements(tempt,"#curve_z#",z_f.replace("_","\\string_"))

    controlForces(log_f)
    keyframeDiff(log_f)
    Ef_add_Ek(log_f)

    energy_f = opt_fig(log_f)
    changeElements(tempt,"#inner_it#",energy_f[0].replace("_","\\string_"))
    changeElements(tempt,"#outer_it#",energy_f[1].replace("_","\\string_"))

    changeElements(tempt,"#number_of_nodes#", grepInt(log_f, "obj vertex number"))
    changeElements(tempt,"#number_of_tets#", grepInt(log_f, "vol tet number"))
    changeElements(tempt,"#T#", grepInt(log_f, "total frame number") )
    changeElements(tempt,"#time_step#", grepFloat(log_f, "h:") )
    changeElements(tempt,"#rs#", grepFloat(log_f, "rs:") )
    changeElements(tempt,"#rw#", grepFloat(log_f, "rw:") )
    changeElements(tempt,"#alphak#", grepFloat(log_f, "alpha_k:") )
    changeElements(tempt,"#alpham#", grepFloat(log_f, "alpha_m:") )
    changeElements(tempt,"#con_nodes#", grepFloat(log_f, "num_of_con_nodes:") )
    changeElements(tempt,"#con_frames#", grepFloat(log_f, "num_of_con_frames:") )
    changeElements(tempt,"#lambda0#", grepFloatList(log_f, "lambda0") )
    changeElements(tempt,"#lambda#", grepFloatList(log_f, "new_lambda") )
    changeElements(tempt,"#optimal fun#", grepFloat(log_f, "final optimal value:") )
    changeElements(tempt,"#opt_time#", grepFloat(log_f, "interpolate time:") )
        
    changeElements(tempt,"#number_of_outer_its#",str(count(log_f,"outer-it-fun")))
    changeElements(tempt,"#number_of_it(Z)#",str(count(log_f,"inner-it-fun(Z)")))
    changeElements(tempt,"#number_of_it(S)#",str(count(log_f,"inner-it-fun-grad(S)")))
    changeElements(tempt,"#number_of_it(l)#",str(count(log_f,"inner-it-fun(Lambda)")))
    changeElements(tempt,"#number_of_it(K)#",str(sumInt(log_f,"Number of Iterations....:")))

    changeElements(tempt,"#lbfgs#",sum(log_f,"lbfgs"))
    changeElements(tempt,"#lbfgs_val_grad#",sum(log_f,"FunGradOptS"))

    changeElements(tempt,"#ipopt_K#",sum(log_f,"MtlOptSolver::optimizeAtA_ipopt()"))

    changeElements(tempt,"#sqp#",sum(log_f,"optimizeZ"))
    changeElements(tempt,"#prepare_for_new_x#",sum(log_f,"prepare_for_new_x"))
    changeElements(tempt,"#sqp_val#",sum(log_f,"MtlOptSqpFunction::val"))
    changeElements(tempt,"#sqp_grad#",sum(log_f,"MtlOptSqpFunction::gra"))
    changeElements(tempt,"#sqp_hess#",sum(log_f,"MtlOptSqpFunction::hes"))
    changeElements(tempt,"#linear_solver::create#",sum(log_f,"linear_solver::create"))
    changeElements(tempt,"#linear_solver::solve#",sum(log_f,"linear_solver::solve"))
    changeElements(tempt,"#EcZ_jacbian#",sum(log_f,"PosConEnergyZ::jacbian"))

    changeElements(tempt,"#sqp_NP#",str(count(log_f,"NP point")))
    changeElements(tempt,"#sqp_TC#",str(count(log_f,"TC point")))    
    changeElements(tempt,"#sqp_DL#",str(count(log_f,"DL point")))

    initfilename = grepStr(log_f,"begin to load initfile:")
    if len(initfilename) > 0:
        changeElements(tempt,"#init_file_name#", initfilename.split()[-1].replace("_","\\_") )
    else:
        changeElements(tempt,"#init_file_name#", "" )

    log_f = "./tempt_opt/"+log_f
    tex_str += open(tempt).read()

open(report_tex,'w').write(tex_str+"\n\end{document}")
os.system("pdflatex "+" -output-directory="+os.path.dirname(report_tex)+" "+report_tex)
