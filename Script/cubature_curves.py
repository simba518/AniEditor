#! /usr/bin/env python
import matplotlib.pyplot as plt

cub = open("./Script/cubature_curves")

error = [1.0]
samples = [0]

for line in cub:
    samples.append(int(line.split(" ")[2]))
    error.append(float(line.split(" ")[6].split("\n")[0].split("%")[0])/100.0)

plt.plot(samples,error,"-ro",label="cubatures",linewidth=5,markersize=10)
plt.ylabel('Relative error')
plt.xlabel('Cubature number')
plt.savefig("./Script/cubature_curves.png")
