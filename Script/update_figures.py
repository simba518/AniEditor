#! /usr/bin/env python

import os
from resize_images import *
from utility import *

paper_doc = "/home/simba/Workspace/AnimationEditorDoc/paper/"
image_doc = paper_doc+"/images/"
model = "beam_fine"

# control forces
Ef = grepFloatNumbers("./tempt_opt/"+model+"-mtlopt.ini.mtllog-Ef-mtlopt.png.txt")
Ef_no = grepFloatNumbers("./tempt_opt/"+model+"-mtlopt.ini.mtllog-Ef.png.txt")
plt.plot(range(0,len(Ef)/2),Ef[0:len(Ef)/2],linewidth=5,label='Material optimization')
plt.plot(range(0,len(Ef_no)/2),Ef_no[0:len(Ef)/2],linewidth=5,label='No material optimization')
legend = plt.legend()
for label in legend.get_texts():
    label.set_fontsize(30)
plt.xlabel("Frame number",fontsize=30)
plt.ylabel("$E_f$",fontsize=30)
fig = plt.gcf()
fig.set_size_inches(22,8.5)
plt.savefig(image_doc+"/control_froces.png")

# curve z
model = "beam_fine"
os.system("cp "+"./tempt_opt/"+model+"-mtlopt-ini-mtllog-z.png "+image_doc+"beam_z_no_mtlopt.png")
os.system("cp "+"./tempt_opt/"+model+"-mtlopt-ini-mtllog-z_mtlopt.png "+image_doc+"beam_z.png")

# compile
os.chdir(paper_doc)
os.system("pdflatex "+"./mat_opt.tex")
