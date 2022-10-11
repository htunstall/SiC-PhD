# -*- coding: utf-8 -*-
import os
import rdf
import sys
import glob
import ase.io
import argparse
import pandas            as pd
import seaborn           as sns
import numpy             as np
import networkx          as nx
import matplotlib.pyplot as plt
from ase.neighborlist    import NeighborList
from matplotlib.ticker   import MaxNLocator
from matplotlib.ticker   import MultipleLocator
from itertools           import cycle
from tqdm                import tqdm

print("ASE VERSION:", ase.__version__)
#------------------------------------------------------------------------------

# Bonding
# C - C : 1.6 Angstrom
# Si- C : 1.89 Angstrom
# Si-Si : 2.35 Angstrom

parser = argparse.ArgumentParser()
parser.add_argument("type",                                               type=str, help="The type of directory you wish to plot e.g. 3000K, quench_010ps etc.")
parser.add_argument("dirs",               nargs="+",                      type=str, help="The directories you wish to plot the ring histogram for")
parser.add_argument("--dft",                         action="store_true",           help="Plot DFT data")
parser.add_argument("--elements",         nargs="+", default=["C","Si"],  type=str, help="The elements we wish to plot for comparison")
parser.add_argument("--max_ring_display", nargs="?", default=10,          type=int, help="The maximum size of ring that will be shown on the graph")

args = parser.parse_args()

if args.type.startswith("quench"):
    filename = "dump.300K-equil.nc"
elif args.type == "3000K":
    filename = "dump.3000K-equil.nc"
else:
    sys.exit("Filename cannot be infered from the type of test")

dft_dir = None
if args.dft:
    #Plot DFT first
    if args.type.startswith("quench"):
        dft_dir      = "/storage/chem/msufgx/postgrad/software/SiC-framework/testing-dir/CP2K/long-quench_10ps"
        dft_filename = "300K-step_mod.SiC-pos-1.xyz"
    elif args.type == "3000K":
        dft_dir = "/storage/chem/msufgx/postgrad/software/SiC-framework/testing-dir/CP2K/long_NPT_10ps"
        dft_filename = "mod.SiC-pos-1.xyz"
    else:
        sys.exit("We don't have dft data for that tpye of simulation")

if dft_dir is not None:
    args.dirs.insert(0, dft_dir)


rings_per_dir = []
for i, _dir in enumerate(args.dirs):
    old_root = os.getcwd()

    if i == 0 and dft_dir is not None:
        old_filename = filename
        filename = dft_filename

        os.chdir(_dir)
    else:
        os.chdir(os.path.join(_dir, args.type))


    all_rings_dict = {}
    for element in args.elements:
        data_filename = "rings_stride-*_element-{}.txt".format(element)
        data_filename = glob.glob(data_filename)[0]


        if os.path.isfile(data_filename):
            all_rings_dict[element] = np.loadtxt(data_filename)
    
        else:
            sys.exit("Can't read ring data in directory: {}".format(_dir))
    os.chdir(old_root)

    rings_per_dir.append(all_rings_dict)

    if i == 0 and dft_dir is not None:
        filename = old_filename

#------------------------------------------------------------------------------
# Plotting
#------------------------------------------------------------------------------
fontsize    = 16

colours     = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
styles      = ["-", "-.", "--", ":"]
linecyclera = cycle(["k"] + [c+s for s in styles for c in colours])
linecyclerb = cycle([c+s for s in styles for c in colours])


#------------------------------------------------------------------------------
# Rings Hist
#width     = 0.7 / len(args.dirs)
fig, axes = plt.subplots(2,1, figsize=(20,15))

for plot_i, element in enumerate(args.elements):
    df_data = []
    for i, (all_rings_dict, _dir) in enumerate(zip(rings_per_dir, args.dirs)):
        data = [(index, value, _dir, element) for index, _list in enumerate(all_rings_dict[element].T) if index < args.max_ring_display for value in _list]
        for item in data:
            df_data.append(item)



    df = pd.DataFrame(df_data, columns=["index", "Count", "Location", "Species"])


    ax = sns.barplot(x="index", y="Count", hue="Location", data=df, ax=axes[plot_i], estimator=np.mean, ci=68, capsize=0.2) # #ci="sd"



axes[0].set_title("Average Number of Rings per Rings Size (Carbon)",  fontsize=fontsize+2)
axes[1].set_title("Average Number of Rings per Rings Size (Silicon)", fontsize=fontsize+2)

for ax in axes:
    ax.set_xlabel("Rings Size", fontsize=fontsize)
    ax.set_ylabel("Average Number of Rings", fontsize=fontsize)
    
    ax.tick_params(axis="both", labelsize=fontsize)
    ax.set_yscale("log")
#    ax.legend(fontsize=fontsize)
    ax.set_axisbelow(True)
    ax.grid()
    
#    ax.xaxis.set_minor_locator(MultipleLocator(1))
#    ax.set_xlim(0, max_chain)
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=2.5)

fig.tight_layout(pad=3)

fig.savefig("combined-bar-plot_{}_ring-analysis_cut-{}.png".format(args.type, args.max_ring_display), dpi=300, bbox_inches="tight")
