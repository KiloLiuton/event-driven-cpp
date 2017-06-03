import os
import numpy as np
import matplotlib.pyplot as plt

dataDir = "relaxationData/"
allfiles = os.listdir(dataDir)
datafiles = sorted([f for f in allfiles if ".txt" in f])
for i,f in enumerate(datafiles): print(i,f)
indexes = [int(i) for i in input("select files for opening:").split(" ")]
smoothdata = int(input("enter amount of smoothing (0 for none):"));

for idx in indexes:
    data = np.loadtxt(dataDir+datafiles[idx], skiprows=2)

    dt = data[:,0]
    dtsum = 0
    for i in range(len(dt)):
        dtsum += dt[i]
        dt[i] = dtsum

    r = data[:,1]

    if(smoothdata):
        newdt = np.empty(shape=(0,))
        newr = np.empty(shape=(0,))
        newlen = int(len(dt) / smoothdata)
        for i in range(newlen):
            j = smoothdata * i
            dtslice = dt[j : j + smoothdata]
            rslice = r[j : j + smoothdata]
            newdt = np.append(newdt, np.average(dtslice))
            newr = np.append(newr, np.average(rslice))
    else:
        newdt = dt
        newr = r
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    ax.plot(newdt, newr, 'r-', lw=0.8)
plt.show()
