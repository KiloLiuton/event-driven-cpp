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

ax.set_title("|r| vs a")
ax2.set_title("X vs a")

for i,idx in enumerate(indexes):
    f = datafiles[idx]
    data = np.loadtxt(dataDir+f, skiprows=2)

    a = data[:,0]
    r = data[:,1]
    X = data[:,2]
    Xnew = data[:,3]

    Xmax = np.max(X)
    Xnewmax = np.max(Xnew)

    X = X / Xmax
    Xnew = Xnew / Xnewmax

    # get N, k from filename
    lab = "N=" + f[f.find("N=")+2 : f.find("k=")] + " k=" + f[f.find("k=")+2 : f.find("p=")]
    print(lab)
    
    ax.plot(a, r,'o', color=plot_colors[i], label=lab)
    ax2.plot(a, Xnew, '-', color=plot_colors[i], lw=1)
    ax.legend(numpoints=1, loc='best')
plt.show()
