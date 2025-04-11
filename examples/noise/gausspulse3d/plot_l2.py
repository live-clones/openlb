import matplotlib.pyplot as plt
import numpy as np
from glob import glob

plt.rcParams['text.usetex'] = True

fig, ax = plt.subplots(figsize=(6,3), dpi=300)
filename = glob("tmp_eternal_*_res101*/gnuplotData/data/l2_absolute.dat")[0]
data = np.loadtxt(filename)
ax.plot(data[:,0], data[:,1], label="Idealfall (vergrößerte Domäne)")
filename = glob("tmp_periodic_*_res101*/gnuplotData/data/l2_absolute.dat")[0]
data = np.loadtxt(filename)
ax.plot(data[:,0], data[:,1], label="Periodische Randbedingungen")
filename = glob("tmp_local_*_res101*/gnuplotData/data/l2_absolute.dat")[0]
data = np.loadtxt(filename)
ax.plot(data[:,0], data[:,1], label="Lokale Randbedingungen")
for bd in [5, 10, 20, 30, 50]:
    filename = glob("tmp_damping_*_res101_bd" + str(bd) + "*/gnuplotData/data/l2_absolute.dat")[0]
    data = np.loadtxt(filename)
    ax.plot(data[:,0], data[:,1], label=f"Dämpfungsschicht {bd} Punkte")
ax.set_yscale('log')
ax.set_ylim((1e-3,1.1))
ax.set_ylabel(r"$\frac{L_p}{L_{p0}}$", rotation=0, ha='right', fontsize=15)
ax.set_xlabel(r"$t_\mathrm{LU}$", fontsize=11, va='top')
ax.legend()
plt.tight_layout()
plt.show()
fig.savefig("l2ComparisonPlot.png")