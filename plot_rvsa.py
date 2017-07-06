import os
import numpy as np
import matplotlib.pyplot as plt
import colorsys

def getColorList(N):
    HSV_colors = [(i/N, 0.5, 1.0) for i in range(N)]
    RGB_colors = [colorsys.hsv_to_rgb(hsv_color[0], hsv_color[1], hsv_color[2]) for hsv_color in HSV_colors]
    return RGB_colors

dataDir = "rvsaData/"
allfiles = os.listdir(dataDir)
datafiles = sorted([f for f in allfiles if ".txt" in f])
for i,f in enumerate(datafiles): print(i,f)
indexes = [int(i) for i in input("select files for opening:").split(" ")]
plot_colors = getColorList(len(indexes))

fig = plt.figure()
ax = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
ax2.set_yscale('log')

for i,idx in enumerate(indexes):
    data = np.loadtxt(dataDir+datafiles[idx], skiprows=2)

    a = data[:,0]
    r = data[:,1]
    X = data[:,2]
    Xnew = data[:,3]

    Xmax = np.max(X)
    Xnewmax = np.max(Xnew)

    X = X / Xmax
    Xnew = Xnew / Xnewmax
    
    ax.plot(a, r,'o', color=plot_colors[i])
    ax2.plot(a, Xnew, '-', color=plot_colors[i], lw=1)
plt.show()
