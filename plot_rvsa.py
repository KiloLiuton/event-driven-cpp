import os
import numpy as np
import matplotlib.pyplot as plt

dataDir = "rvsaData/"
allfiles = os.listdir(dataDir)
datafiles = sorted([f for f in allfiles if ".txt" in f])
for i,f in enumerate(datafiles): print(i,f)
indexes = [int(i) for i in input("select files for opening:").split(" ")]

for idx in indexes:
    data = np.loadtxt(dataDir+datafiles[idx], skiprows=2)

    a = data[:,0]
    r = data[:,1]
    X = data[:,2]
    Xnew = data[:,3]

    Xmax = np.max(X)
    Xnewmax = np.max(Xnew)

    X = X / Xmax
    Xnew = Xnew / Xnewmax

    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    
    ax.plot(a, r,'ro')
    ax2.plot(a, X, 'b-', lw=1)
    ax2.plot(a, Xnew, 'r-', lw=1)
plt.show()
