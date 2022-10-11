import rdf
import os
import ase
import sys
import math
import argparse
import matplotlib.pyplot as plt
import numpy             as np

#=================================================================
def plot_rdf(midpoints, data, label, axes):
    """
    data is a list of y values for the plotting
    """
    # Plot in 2d grid
    for i, y in enumerate(data):
        row = math.floor(i / 2)
        col = i % 2
        axes[row,col].plot(midpoints, y, label=label)

#=================================================================
fontsize = 26
verbose  = True



parser = argparse.ArgumentParser()
parser.add_argument("filename",                                   type=str, help="The NetCDF trajectory file")
parser.add_argument("type",                                       type=str, help="The type of directory you wish to plot e.g. 3000K, quench_010ps etc.")
parser.add_argument("dirs",             nargs="+",                type=str, help="The directories you wish to plot the RDF for")
parser.add_argument("--labels",         nargs="+", default=None,  type=str, help="The labels for the RDF plots")
parser.add_argument("-s", "--stride",   nargs="?", default=100,   type=int, help="The stride to take the trajectory at")
parser.add_argument("-d", "--dft",      nargs="?", default="npt", type=str, help="Which DFT are we plotting: NPT, NVT or Amorphous (case insensitive)")

parser.add_argument("-e", "--exp",      action="store_true",                help="Plot the experemntal RDFs")

args = parser.parse_args()

if args.labels is None:
    args.labels = [None] * len(args.dirs)
elif type(args.labels) is list:
    if len(args.labels) != len(args.dirs):
        sys.exit("You have a diffrent number of labels ({}) and directories ({})".format(len(args.labels), len(args.dirs)))

root = os.getcwd()

dfts = ["npt", "nvt", "amorphous"]
if args.dft.lower() not in dfts:
    sys.exit("You have set an incorrect DFT parameter") 

print("Plotting with:\n DFT: {}\n Experemental: {}\n Simulation Type: {}\n".format(args.dft, args.exp, args.type))

# The figures
fig1, axes1 = plt.subplots(2,2, figsize=(20,20)) # total average
   
for _dir, label in zip(args.dirs, args.labels):
    print("Working on: `{}`".format(_dir))
    os.chdir(os.path.join(_dir, args.type))
    # Only lord the dir if we don't have rdf files
    if not os.path.isfile("midpoints_species-all_bins-200_{}.txt".format(args.type)):
        print("Loading the trajectory")
        atoms_list = rdf.lammps2atoms(ase.io.NetCDFTrajectory(args.filename, "r")[::args.stride])
        print("The size is:", len(atoms_list))
    else:
        atoms_list = []

    comment  = os.path.basename(os.getcwd())
    #==============================================
    # Average
    #==============================================
    #                                             cutoff, s1, s2
    cc_array_2d,   midpoints = rdf.get_rdf(atoms_list, 6,  6,  6, verbose=verbose, save=True, comment=comment)
    sisi_array_2d, midpoints = rdf.get_rdf(atoms_list, 6, 14, 14, verbose=verbose, save=True, comment=comment)
    sic_array_2d,  midpoints = rdf.get_rdf(atoms_list, 6, 14,  6, verbose=verbose, save=True, comment=comment)
    all_array_2d,  midpoints = rdf.get_rdf(atoms_list, 6,         verbose=verbose, save=True, comment=comment)
    
    data = [np.mean(cc_array_2d,   axis=0),
            np.mean(sisi_array_2d, axis=0),
            np.mean(sic_array_2d,  axis=0),
            np.mean(all_array_2d,  axis=0)]
    
    if label is not None:
        comment = label
    else:
        comment = _dir

    plot_rdf(midpoints, data, comment, axes1)

    os.chdir(root)
    print()


#==================================================
# NEW DFT
#==================================================
print("Working on: `DFT {}`".format(args.dft))
#--------------------------------------------------
# NPT
#/storage/chem/msufgx/postgrad/software/SiC-framework/testing-dir/lammps/rdf-size-time-effect/structures/216_atoms/dft-like-quench/DFT_NPT
#--------------------------------------------------
if args.dft.lower() == "npt":
    old_root = os.getcwd()
    label = "DFT @ 3000 K (NPT)"
    os.chdir("/storage/chem/msufgx/postgrad/software/SiC-framework/testing-dir/lammps/rdf-size-time-effect/structures/216_atoms/dft-like-quench/DFT_NPT")
    dft_cc_array_2d,   dft_midpoints = rdf.get_rdf([], 6,  6,  6, verbose=verbose, save=True, comment="DFT_NPT")
    dft_sisi_array_2d, dft_midpoints = rdf.get_rdf([], 6, 14, 14, verbose=verbose, save=True, comment="DFT_NPT")
    dft_sic_array_2d,  dft_midpoints = rdf.get_rdf([], 6, 14,  6, verbose=verbose, save=True, comment="DFT_NPT")
    dft_all_array_2d,  dft_midpoints = rdf.get_rdf([], 6,         verbose=verbose, save=True, comment="DFT_NPT")
    os.chdir(old_root)

#--------------------------------------------------
# NVT
#/storage/chem/msufgx/postgrad/software/SiC-framework/testing-dir/lammps/rdf-size-time-effect/structures/216_atoms/dft-like-quench/DFT_NPT
#--------------------------------------------------
elif args.dft.lower() == "nvt":
    old_root = os.getcwd()
    label = "DFT @ 3000 K (NVT)"
    os.chdir("/storage/chem/msufgx/postgrad/software/SiC-framework/testing-dir/lammps/rdf-size-time-effect/structures/216_atoms/dft-like-quench/DFT_NVT")
    dft_cc_array_2d,   dft_midpoints = rdf.get_rdf([], 6,  6,  6, verbose=verbose, save=True, comment="DFT_NVT")
    dft_sisi_array_2d, dft_midpoints = rdf.get_rdf([], 6, 14, 14, verbose=verbose, save=True, comment="DFT_NVT")
    dft_sic_array_2d,  dft_midpoints = rdf.get_rdf([], 6, 14,  6, verbose=verbose, save=True, comment="DFT_NVT")
    dft_all_array_2d,  dft_midpoints = rdf.get_rdf([], 6,         verbose=verbose, save=True, comment="DFT_NVT")
    os.chdir(old_root)

#--------------------------------------------------
# DFT - amorphous
#--------------------------------------------------
elif args.dft.lower() == "amorphous":
    old_root = os.getcwd()
    label = "DFT @ 300 K"
    os.chdir("/storage/chem/msufgx/postgrad/software/SiC-framework/testing-dir/CP2K/long-quench-last-step-300K")
    dft_cc_array_2d,   dft_midpoints = rdf.get_rdf([], 6,  6,  6, verbose=verbose, save=True, comment="long-quench-last-step-300K")
    dft_sisi_array_2d, dft_midpoints = rdf.get_rdf([], 6, 14, 14, verbose=verbose, save=True, comment="long-quench-last-step-300K")
    dft_sic_array_2d,  dft_midpoints = rdf.get_rdf([], 6, 14,  6, verbose=verbose, save=True, comment="long-quench-last-step-300K")
    dft_all_array_2d,  dft_midpoints = rdf.get_rdf([], 6,         verbose=verbose, save=True, comment="long-quench-last-step-300K")
    os.chdir(old_root)


#        location, y-axis data,                        label, title bonding
to_plot = [[(0,0), np.mean(dft_cc_array_2d,   axis=0), label, "C-C"  ],
           [(0,1), np.mean(dft_sisi_array_2d, axis=0), label, "Si-Si"],
           [(1,0), np.mean(dft_sic_array_2d,  axis=0), label, "Si-C" ],
           [(1,1), np.mean(dft_all_array_2d,  axis=0), label, "All"  ]]

#==================================================
# Plot in 2d grid
#==================================================
for (i,j), y, label, title in to_plot:
    axes1[i,j].plot(midpoints, y, label=label, color="k")


#==================
# Experemental
#==================
if args.exp:
    exp_data = np.loadtxt("/storage/chem/msufgx/postgrad/software/SiC-framework/analysis/RDF/sorted.experemental_800C-annealed-amorphous.txt")
    axes1[1,1].plot(exp_data[:,0]*10, exp_data[:,1], "--", c="k", label="annealed 800C (Exp.)")
    
    exp_data = np.loadtxt("/storage/chem/msufgx/postgrad/software/SiC-framework/analysis/RDF/filtered.csv", delimiter=",")
    axes1[1,1].plot(exp_data[:,0]*10, exp_data[:,1], "-.", c="k", label="Filtered (Exp.)")
    
    exp_data = np.loadtxt("/storage/chem/msufgx/postgrad/software/SiC-framework/analysis/RDF/unfiltered.csv", delimiter=",")
    axes1[1,1].plot(exp_data[:,0]*10, exp_data[:,1], ".", c="k", label="Unfiltered (Exp.)")


#==================================================
# Formatting
#==================================================
titles = [[(0,0), "C-C"  ],
          [(0,1), "Si-Si"],
          [(1,0), "Si-C" ],
          [(1,1), "All"  ]]

for (i,j), title in titles:
    axes1[i,j].set_title("RDF for {} bonds".format(title), fontsize=fontsize+2)

for ax in axes1.flatten():
    ax.tick_params(axis="both", which="major", labelsize=fontsize-2)
    ax.set_xlabel("Distance [$\AA$]", fontsize=fontsize)
    ax.set_ylabel("g(r)",             fontsize=fontsize)
#    ax.legend(ncol=3, bbox_to_anchor=(0.5, -0.1), loc="upper center", fontsize=fontsize-4)
    ax.grid()

# Set global legend
handles, labels = axes1[1,1].get_legend_handles_labels()
lgd = fig1.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.01), fontsize=fontsize-4, ncol=3)


fig1.tight_layout(pad=3)
fig1.subplots_adjust(bottom=0.2)
#==================================================
# Saving
#==================================================
fig1.savefig("rdf_{}_combined.png".format(args.type), dpi=150, bbox_inches="tight", bbox_extra_artists=(lgd,))
fig1.savefig("rdf_{}_combined_publication.png".format(args.type), dpi=600, bbox_inches="tight", bbox_extra_artists=(lgd,), transparent=True)

print("Finished!")
