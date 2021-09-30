import os
import sys
import rdf
import argparse
import ase.io
import numpy               as np
import matplotlib.pyplot   as plt
from   matscipy.rings  import ring_statistics
from matplotlib.ticker import MaxNLocator
from tqdm              import tqdm


parser = argparse.ArgumentParser()
parser.add_argument("filename", metavar="f", type=str, help="The NetCDF trajectory file")
parser.add_argument("stride",   metavar="s", nargs="?", const=100, type=int, help="The stride to take the trajectory at")

args = parser.parse_args()



atoms_list = rdf.lammps2atoms(ase.io.NetCDFTrajectory(args.filename, "r")[::args.stride])
#    atoms_list   = ase.io.read("rings.xyz", index=":")


if not os.path.isfile("rings_stride{}.txt".format(args.stride)):

    hist_2d = np.zeros((len(atoms_list),25))
    for i, atoms in enumerate(tqdm(atoms_list)):
        # Remove the silicon atoms
        del atoms[[atom.index for atom in atoms if atom.symbol=="Si"]]

        data = ring_statistics(atoms, cutoff=2.0)
        for j in range(len(data)):
            hist_2d[i,j] = data[j]

    print()
    np.savetxt("rings_stride{}.txt".format(args.stride), hist_2d)


else:
    hist_2d = np.loadtxt("rings_stride{}.txt".format(args.stride))


fontsize = 16
fig, ax  = plt.subplots(figsize=(20,10))

time = np.linspace(float(atoms_list[0].info["time"]), float(atoms_list[-1].info["time"]), len(atoms_list))

for i in range(3, 11, 1):
    if hist_2d[:,i].sum() != 0:
        ax.plot(time, hist_2d[:,i], label= "{}".format(i))
    else:
        print("No {} member rings, not plotting".format(i))


##------------------------------------------------------------------
##colvar = np.loadtxt("../COLVAR") # time, cn.mean, meta.bias, uwall.bias
#data = np.loadtxt("std.out", skiprows=26, max_rows=20001) # temp = data[:,2]
#
#ax.yaxis.set_major_locator(MaxNLocator(integer=True))
#ax2 = ax.twinx()
##colvar[:,0] = colvar[:,0] - colvar[0,0]
#
#ax2.plot(time,data[:,2], color="k")
#
#ax2.set_ylabel("Temperature [K]", fontsize=fontsize)
#
##ax2.set_ylim(ax2.get_ylim()[0]-ax2.get_ylim()[0]*.2, ax2.get_ylim()[1])
##------------------------------------------------------------------





ax.legend(fontsize=fontsize, fancybox=True, shadow=True, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.1), title="Ring Size") #, title_fontsize=16)
ax.tick_params(axis="both", labelsize=16)

ax.yaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_ylabel("Number of rings", fontsize=fontsize)
ax.set_xlabel("Time [ps]", fontsize=fontsize)
ax.grid()


ax.set_xlim((time[0], time[-1]))
#ax.set_ylim((0, ax.get_ylim()[1]+ax.get_ylim()[1]*0.2))

#ax.set_xlim(21.5,24)

fig.tight_layout(pad=15)
fig.savefig("rings_stride{}.png".format(args.stride), dpi=300)
