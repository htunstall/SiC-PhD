# -*- coding: utf-8 -*-
import os
import rdf
import sys
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
parser.add_argument("type",                                       type=str, help="The type of directory you wish to plot e.g. 3000K, quench_010ps etc.")
parser.add_argument("dirs",             nargs="+",                type=str, help="The directories you wish to plot the chain histogram for")
parser.add_argument("-d", "--dft",      action="store_true",                help="Plot DFT data")

args = parser.parse_args()

species         = ["C", "Si"]
graph_max_chain = 15

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


max_chain = 45
chains_per_dir = []
for i, _dir in enumerate(args.dirs):
    old_root = os.getcwd()

    if i == 0 and dft_dir is not None:
        old_filename = filename
        filename = dft_filename

        os.chdir(_dir)
    else:
        os.chdir(os.path.join(_dir, args.type))


    all_chains_dict = {}
    for element in species:
        if os.path.isfile("{}_chains_{}.txt".format(element.lower(), ".".join(filename.split(".")[:-1]))):
            all_chains_dict[element] = np.loadtxt("{}_chains_{}.txt".format(element.lower(), ".".join(filename.split(".")[:-1])))
    
        else:
            sys.exit("Can't read chain data in directory: {}".format(_dir))#
    os.chdir(old_root)

    chains_per_dir.append(all_chains_dict)

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
# Chains Hist
#width     = 0.7 / len(args.dirs)
fig, axes = plt.subplots(2,1, figsize=(20,15))

for plot_i, species in enumerate(species):
    df_data = []
    for i, (all_chains_dict, _dir) in enumerate(zip(chains_per_dir, args.dirs)):
        data = [(index, value, _dir, species) for index, _list in enumerate(all_chains_dict[species].T) if index < graph_max_chain for value in _list]
        for item in data:
            df_data.append(item)



    df = pd.DataFrame(df_data, columns=["index", "Count", "Location", "Species"])


    ax = sns.barplot(x="index", y="Count", hue="Location", data=df, ax=axes[plot_i], estimator=np.mean, ci=68, capsize=0.2) # #ci="sd"



#    x        = np.arange(1, max_chain, 1, dtype=int)
#    c_means  = np.mean(all_chains_dict["C"], axis=0)
#    c_stdds  = np.std(all_chains_dict["C"], axis=0)
#    
#    si_means = np.mean(all_chains_dict["Si"], axis=0)
#    si_stdds = np.std(all_chains_dict["Si"], axis=0)
#    
#    
#    c_bars  = ax[0].bar(x - width/2, c_means[1:],  width, label="Carbon {}".format(_dir))
#    si_bars = ax[1].bar(x + width/2, si_means[1:], width, label="Silicon {}".format(_dir))
#    
#    # Error Bars
#    ax[0].errorbar(x - width/2, c_means[1:],  yerr=c_stdds[1:],  ecolor="k", barsabove=True, fmt="none")
#    ax[1].errorbar(x + width/2, si_means[1:], yerr=si_stdds[1:], ecolor="k", barsabove=True, fmt="none")



axes[0].set_title("Average Number of Chains per Chain Length (Carbon)", fontsize=fontsize+2)
axes[1].set_title("Average Number of Chains per Chain Length (Silicon)", fontsize=fontsize+2)

for ax in axes:
    ax.set_xlabel("Chain Length", fontsize=fontsize)
    ax.set_ylabel("Average Number of Chains", fontsize=fontsize)
    
    ax.tick_params(axis="both", labelsize=fontsize)
    ax.set_yscale("log")
#    ax.legend(fontsize=fontsize)
    ax.set_axisbelow(True)
    ax.grid()
    
#    ax.xaxis.set_minor_locator(MultipleLocator(1))
#    ax.set_xlim(0, max_chain)
    ax.set_ylim(bottom=0)

fig.tight_layout(pad=3)

fig.savefig("combined-bar-plot_chain-analysis_cut-{}.png".format(graph_max_chain), dpi=300, bbox_inches="tight")
