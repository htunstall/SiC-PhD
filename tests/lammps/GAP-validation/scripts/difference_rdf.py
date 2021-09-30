import rdf
import os
import ase
import sys
import math
import argparse
import matplotlib.pyplot as plt
import numpy             as np

def load_rdf(comment, root, traj_name, args):
    old_root = os.getcwd()
    os.chdir(root)

    # If the RDF has been run
    if os.path.isfile("array-2d_species-all_bins-200_{}.txt".format(comment)):
        al = []
    else:
        al = rdf.lammps2atoms(ase.io.NetCDFTrajectory(traj_name, "r")[::args.stride])

    cc_array_2d,   midpoints = rdf.get_rdf(al, 6,  6,  6, verbose=verbose, save=True, comment=comment)
    sisi_array_2d, midpoints = rdf.get_rdf(al, 6, 14, 14, verbose=verbose, save=True, comment=comment)
    sic_array_2d,  midpoints = rdf.get_rdf(al, 6, 14,  6, verbose=verbose, save=True, comment=comment)
    all_array_2d,  midpoints = rdf.get_rdf(al, 6,         verbose=verbose, save=True, comment=comment)
    os.chdir(old_root)
    
    values = [np.mean(cc_array_2d,   axis=0),
              np.mean(sisi_array_2d, axis=0),
              np.mean(sic_array_2d,  axis=0),
              np.mean(all_array_2d,  axis=0)]

    return values, midpoints

#---------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("type",                                         type=str, help="The type of directory you wish to plot e.g. 3000K, quench_010ps etc.")
parser.add_argument("dirs",               nargs="+",                type=str, help="The directories you wish to plot the RDF for")
parser.add_argument("-c", "--comparison", nargs="?", default="npt", type=str, help="What are we comparing to? npt, nvt, amorphous, exp, expu, expf")
parser.add_argument("-s", "--stride",     nargs="?", default=100,   type=int, help="The stride to take the trajectory at")

args = parser.parse_args()
#---------------------------------------------------
verbose  = True
fontsize = 20

types = ["3000K", "quench_010ps"]
if args.type not in types:
    sys.exit("You have set an incorrect simulation type")

comparisons = {"npt"       : "/storage/chem/msufgx/postgrad/software/SiC-framework/testing-dir/lammps/rdf-size-time-effect/structures/216_atoms/dft-like-quench/DFT_NPT",
               "nvt"       : "/storage/chem/msufgx/postgrad/software/SiC-framework/testing-dir/lammps/rdf-size-time-effect/structures/216_atoms/dft-like-quench/DFT_NVT",
               "amorphous" : "/storage/chem/msufgx/postgrad/software/SiC-framework/testing-dir/CP2K/long-quench-last-step-300K"}
# NOT YET IMPLEMENTED
#               "exp"       : "/storage/chem/msufgx/postgrad/software/SiC-framework/analysis/RDF/sorted.experemental_800C-annealed-amorphous.txt",
#               "expu"      : "/storage/chem/msufgx/postgrad/software/SiC-framework/analysis/RDF/unfiltered.csv",
#               "expf"      : "/storage/chem/msufgx/postgrad/software/SiC-framework/analysis/RDF/filtered.csv"}

labels_dict = {"npt"       : "DFT @ 3000 K (NPT)",
               "nvt"       : "DFT @ 3000 K (NVT)",
               "amorphous" : "DFT @ 300 K"}


if args.comparison.lower() not in comparisons.keys():
    sys.exit("You have set an incorrect DFT parameter")



print("Plotting with:\n Comparion to:    {}\n Simulation Type: {}\n".format(args.comparison, args.type))

if args.type.lower() == "3000k":
    traj_name = "dump.3000K-equil.nc"
elif args.type.lower().startswith("quench"):
    traj_name = "dump.300K-equil.nc"

#----------------------------------------
# Get the RDF to compare to
                                                                                                          
comp_values, midpoints = load_rdf(os.path.basename(comparisons[args.comparison]), comparisons[args.comparison], traj_name, args)
#----------------------------------------

diff_squareds = []
label_list    = []
for _dir in args.dirs:
    print()
    values, _ = load_rdf(args.type, os.path.join(_dir,args.type), traj_name, args)

    diff_squareds.append([(value - comp)**2 for value, comp in zip(values, comp_values)])
    label_list.append(_dir)


#==================================================
# Plotting and Integrals
#==================================================
fig, axes = plt.subplots(2,2, figsize=(20,20))

integrals = np.zeros((len(args.dirs), 4))
for k, (diff_squared, label) in enumerate(zip(diff_squareds, label_list)):
    #        location, y-axis data,     label
    to_plot = [[(0,0), diff_squared[0], label],
               [(0,1), diff_squared[1], label],
               [(1,0), diff_squared[2], label],
               [(1,1), diff_squared[3], label]]
    
    #==================================================
    # Plot in 2d grid
    #==================================================
    for (i,j), y, label in to_plot:
        axes[i,j].plot(midpoints, y, label=label)
        integrals[k,i*2+j] = np.trapz(y, midpoints)

#==================================================
# Formatting
#==================================================
titles = [[(0,0), "C-C"  ],
          [(0,1), "Si-Si"],
          [(1,0), "Si-C" ],
          [(1,1), "All"  ]]

for (i,j), title in titles:
    axes[i,j].set_title("RDF for {} bonds".format(title), fontsize=fontsize+2)

for ax in axes.flatten():
    ax.tick_params(axis="both", which="major", labelsize=fontsize-2)
    ax.set_xlabel("Distance [$\AA$]",   fontsize=fontsize)
    ax.set_ylabel("Difference Squared", fontsize=fontsize)
#    ax.legend(ncol=3, bbox_to_anchor=(0.5, -0.1), loc="upper center", fontsize=fontsize-4)
    ax.grid()

print("\nIntegrals")
print(integrals)
print()

# Set global legend
handles, labels = axes[1,1].get_legend_handles_labels()
lgd = fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.01), fontsize=fontsize-4)

fig.suptitle("Difference Squared RDF vs {} (for test {})".format(labels_dict[args.comparison], args.type), fontsize=fontsize+2, y=1.05)

fig.tight_layout(pad=3)
fig.subplots_adjust(bottom=0.2)
#==================================================
# Saving
#==================================================
fig.savefig("rdf_diff-squared_{}.png".format(args.type), dpi=300, bbox_inches="tight", bbox_extra_artists=(lgd,))

np.savetxt("integrals_{}.txt".format(args.type), integrals)

with open("integral_locations_{}.txt".format(args.type), "w") as f:
    f.writelines([str(os.path.join(_dir, args.type)) + "\n" for _dir in args.dirs])

##==================
## Experemental
##==================
#if args.exp:
#    exp_data = np.loadtxt("/storage/chem/msufgx/postgrad/software/SiC-framework/analysis/RDF/sorted.experemental_800C-annealed-amorphous.txt")
#    axes1[1,1].plot(exp_data[:,0]*10, exp_data[:,1], "--", c="k", label="annealed 800C (Exp.)")
#
#    exp_data = np.loadtxt("/storage/chem/msufgx/postgrad/software/SiC-framework/analysis/RDF/filtered.csv", delimiter=",")
#    axes1[1,1].plot(exp_data[:,0]*10, exp_data[:,1], "-.", c="k", label="Filtered (Exp.)")
#
#    exp_data = np.loadtxt("/storage/chem/msufgx/postgrad/software/SiC-framework/analysis/RDF/unfiltered.csv", delimiter=",")
#    axes1[1,1].plot(exp_data[:,0]*10, exp_data[:,1], ".", c="k", label="Unfiltered (Exp.)")
#
