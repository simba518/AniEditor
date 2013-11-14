#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 26,
        }

f = open("./Script/tempt_beam_mtl_opt50.txt")
recordObjtive = False
iterLabel = "iter "
Y = [6.1209828e+09]
X = [0]
gX = []
gY = [6.1209828e+09]
gX2 = [0]
for line in f:
    if len(line)>len(iterLabel) and iterLabel==line[0:len(iterLabel)]:
        recordObjtive = True
    elif recordObjtive and len(line) > 4 and line[3] <= '9' and line[3] >= '0':
        y = float(line.split()[1])
        if (len(Y)>0 and Y[-1] < y):
            y = Y[-1]
        X.append( len(Y)-1 )
        Y.append( y )
    else:
        if recordObjtive and len(X) > 0:
            gX.append( X[-1] )
            gY.append( Y[-1] )
            gX2.append( len(gY)/2.0-0.5 )
        recordObjtive = False
f.close()
plt.rcParams.update({'font.size': 26})

#################### draw outer iterations.
plt.plot(gX2,gY,'ro-',label='outer iteration')
plt.ylim(min(gY[1:-1])/1.1,max(gY[1:-1])*1.3)
plt.xlim(-1,gX2[-1]*1.1)
#################### draw outer iterations.


#################### draw inner iterations.
# plt.plot(X[0:-1],Y[0:-1],label='inner iteration')
# plt.ylim(-max(Y[0:-1])/10.0,max(Y[0:-1])*1.3)
# plt.xlim(-10,max(X[0:-1])*1.3)
#################### draw inner iterations.

plt.ylabel('Energy', fontdict=font)
plt.xlabel('Iterations', fontdict=font)
plt.legend()
plt.show()
