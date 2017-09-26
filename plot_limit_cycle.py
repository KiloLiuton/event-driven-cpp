import os
import colorsys
import numpy as np
import matplotlib.pyplot as plt

def getColors(N):
    HSV_colors = [((i+1)/N, 0.5, 1.0) for i in range(N)]
    RGB_colors = [colorsys.hsv_to_rgb(c[0], c[1], c[2]) for c in HSV_colors]
    return RGB_colors

dataDir = "relaxationData/"
allfiles = os.listdir(dataDir)
datafiles = sorted([f for f in allfiles if ".txt" in f])
for i,f in enumerate(datafiles): print(i,f)
indexes = [int(i) for i in input("select files for opening:").split(" ")]

colors = getColors(len(indexes))
filename = datafiles[indexes[0]]
N = int(filename[filename.find("N=")+2 : filename.find("k=")])

fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1,1,1)

ax.axvline(0)
ax.axhline(0)

limx = [0, N]
limy = [N, 0]
ax.plot(limx,limy, "b-")

data = np.loadtxt(dataDir+filename, skiprows=2)
p1 = data[:,2][:-1]
p2 = data[:,3][:-1]

ax.plot(p1,p2,"r-")

plt.show()
