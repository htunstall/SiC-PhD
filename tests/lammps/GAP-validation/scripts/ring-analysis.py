import os
import sys
import rdf
import argparse
import ase.io
import pandas              as pd
import numpy               as np
import seaborn             as sns
import matplotlib.pyplot   as plt
from   matscipy.rings  import ring_statistics
from matplotlib.ticker import MaxNLocator
from tqdm              import tqdm


parser = argparse.ArgumentParser()
parser.add_argument("filename",           metavar="f",                                  type=str,   help="The trajectory file")
parser.add_argument("--stride",           metavar="-s", nargs="?", default=100,         type=int,   help="The stride to take the trajectory at")
parser.add_argument("--elements",         metavar="-e", nargs="+", default=["C","Si"],  type=str,   help="The elements we wish to perform ring analysis on")
parser.add_argument("--element_cutoffs",  metavar="-c", nargs="+", default=[1.6, 2.35], type=float, help="The element cutoff values for the provided elements")
parser.add_argument("--max_ring",         metavar="-m", nargs="?", default=25,          type=int,   help="The maximum size of ring that will be stored")
parser.add_argument("--max_ring_display", metavar="-M", nargs="?", default=10,          type=int,   help="The maximum size of ring that will be shown on the graph")

args = parser.parse_args()

atoms_list = None
# Load the Atoms list from either netCDF or extxyz
#  there is no reason why it can't be loaded from any ase reader, just not implemented properly
if args.filename.split(".")[-1] == "nc":
    atoms_list = rdf.lammps2atoms(ase.io.NetCDFTrajectory(args.filename, "r")[::args.stride])
elif args.filename.split(".")[-1] == "xyz":
    atoms_list = ase.io.read(args.filename, index="::{}".format(args.stride))
else:
    sys.exit("Incorrect filename extension, neither `.nc` or `.xyz`")


structure       = [[element, cutoff, [not_element for not_element in args.elements if not_element != element]] for element, cutoff in zip(args.elements, args.element_cutoffs)]

all_rings_dict = {}
for element, cutoff, bad_species in structure:
    if not os.path.isfile("rings_stride-{}_element-{}.txt".format(args.stride, element)):
        all_rings_dict[element] = np.zeros((len(atoms_list),args.max_ring))
        for i, atoms in enumerate(tqdm(atoms_list.copy())):
            # Remove the species we don't want
            del atoms[[atom.index for atom in atoms if atom.symbol in bad_species]]
    
            data = ring_statistics(atoms, cutoff=cutoff)
            for j in range(len(data)):
                all_rings_dict[element][i,j] = data[j]
    
        print()
        np.savetxt("rings_{}_stride-{}_element-{}.txt".format(".".join(args.filename.split(".")[:-1]), args.stride, element), all_rings_dict[element])


    else:
        all_rings_dict[element] = np.loadtxt("rings_stride-{}_element-{}.txt".format(args.stride, element))

#==============================
# Plotting
#==============================
fontsize = 16

#------------------------------------------------------------------------------
# Chains Hist
#------------------------------------------------------------------------------
fig_bar, ax_bar = plt.subplots(1,1, figsize=(10,10))

df_data = []
for element in args.elements:
    data = [(index, value, element) for index, _list in enumerate(all_rings_dict[element].T) if index < args.max_ring_display for value in _list]
    for item in data:
        df_data.append(item)

df     = pd.DataFrame(df_data, columns=["Ring Size", "Count", "Element"])
ax_bar = sns.barplot(x="Ring Size", y="Count", hue="Element", data=df, ax=ax_bar, estimator=np.mean, ci=68, capsize=0.2)


ax_bar.set_title("Average Number of Rings per Ring Size", fontsize=fontsize+2)

ax_bar.set_xlabel("Ring Size", fontsize=fontsize)
ax_bar.set_ylabel("Average Number of Rings", fontsize=fontsize)

ax_bar.tick_params(axis="both", labelsize=fontsize)
ax_bar.set_yscale("log")
ax_bar.set_axisbelow(True)
ax_bar.grid()

ax_bar.set_ylim(bottom=0)

fig_bar.tight_layout(pad=3)

fig_bar.savefig("rings_bar_{}_stride-{}.png".format(".".join(args.filename.split(".")[:-1]), args.stride), dpi=300, bbox_inches="tight")
