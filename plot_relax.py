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
smoothdata = int(input("enter amount of smoothing (0 for none):"));

colors = getColors(len(indexes))

for c,idx in enumerate(indexes):
    data = np.loadtxt(dataDir+datafiles[idx], skiprows=2)

    dt = data[:,0][:-1]
    r = data[:,1][:-1]
    n0 = data[:,2][:-1]
    n1 = data[:,3][:-1]
    relaxationPeriod = int(data[:,0][-1])

    dtsum = 0
    for i in range(len(dt)):
        dtsum += dt[i]
        dt[i] = dtsum

    if(smoothdata):
        newdt = np.empty(shape=(0,))
        newr = np.empty(shape=(0,))
        newn0 = np.empty(shape=(0,))
        newn1 = np.empty(shape=(0,))

        newlen = int(len(dt) / smoothdata)
        for i in range(newlen):
            j = smoothdata * i
            window = dt[j : j + smoothdata]
            newdt = np.append(newdt, np.average(window))

            window = r[j : j + smoothdata]
            newr = np.append(newr, np.average(window))

            window = n0[j : j + smoothdata]
            newn0 = np.append(newn0, np.average(window))

            window = n1[j : j + smoothdata]
            newn1 = np.append(newn1, np.average(window))
    else:
        newdt = dt
        newr = r
        newn0 = n0
        newn1 = n1

    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    ax.plot(newdt, newr, '-', lw=0.8, color=colors[c])
    ax.axvline(dt[relaxationPeriod])
    ax2.plot(newdt, newn0, 'o', color=colors[c])
    ax2.plot(newdt, newn1, 'o', color=colors[(c+1)%len(colors)])
    plt.show()
