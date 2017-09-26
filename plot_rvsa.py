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

fig = plt.figure(figsize=(16,12))
ax = fig.add_axes([0.1, 0.525, 0.6, 0.375])
ax2 = fig.add_axes([0.1, 0.1, 0.6, 0.375])
ax2.set_yscale('log')

FS = 40

f = datafiles[indexes[0]]
N = int(f[f.find("N=")+2 : f.find("k=")])
k = int(f[f.find("k=")+2 : f.find("p=")])
ax.set_title("OP and Variance vs Coup. Str."+" k=" + str(100*(2*k+1)/N)[:5] + "%", fontsize=FS)
ax.set_ylabel("$\langle \langle r \\rangle _t \\rangle _s$", fontsize=FS)
ax2.set_ylabel("$\chi$", fontsize=FS)
ax2.set_xlabel("a", fontsize=FS)

box = ax.get_position()
box2 = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width*0.7, box.height])
#ax2.set_position([box2.x0, box2.y0+box2.height, box2.width*0.7, box2.height])

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

    N = int(f[f.find("N=")+2 : f.find("k=")])
    k = int(f[f.find("k=")+2 : f.find("p=")])
    # get N, k, p from filename
    #lab = "N=" + f[f.find("N=")+2 : f.find("k=")] 
    lab = ""
    #lab += " k=" + str(100*(2*k+1)/N)[:5] + "%"
    lab += " p=" + f[f.find("p=")+2 : f.find("TRIALS=")]
    
    ax.plot(a, r,'-', color=plot_colors[i])
    ax2.plot(a, Xnew, '-', color=plot_colors[i], lw=1, label=lab)
ax2.legend(fontsize=int(0.35*FS), numpoints=1, loc='lower right', bbox_to_anchor=(1.25, 0.6))
plt.savefig("foo.png")
#plt.show()
