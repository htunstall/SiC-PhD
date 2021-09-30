import rdf
import os
import ase
import sys
import math
import argparse
import matplotlib.pyplot as plt
import numpy             as np

#=================================================================
fontsize = 20
verbose  = True
comment  = os.path.basename(os.getcwd())



parser = argparse.ArgumentParser()
parser.add_argument("filename", metavar="f", type=str, help="The NetCDF trajectory file")
parser.add_argument("stride",   metavar="s", nargs="?", default=100,   type=int,  help="The stride to take the trajectory at")
parser.add_argument("dft",      metavar="d", nargs="?", default="npt", type=str,  help="Which DFT are we plotting: NPT, NVT or Amorphous (case insensitive)")
parser.add_argument("exp",      metavar="e", nargs="?", default=False, type=bool, help="Plot the experemntal RDFs")

args = parser.parse_args()


print("Loading the trajectory")
atoms_list = rdf.lammps2atoms(ase.io.NetCDFTrajectory(args.filename, "r")[::args.stride])
print("The size is:", len(atoms_list))
#atoms_list = ase.io.read("mod.SiC-pos-1.xyz", index=":")

dfts = ["npt", "nvt", "amorphous"]
if args.dft.lower() not in dfts:
    sys.exit("You have set an incorrect DFT parameter") 

print("Plotting with:\n DFT: {}\n Experemental: {}".format(args.dft, args.exp))

# The figures
fig1, axes1 = plt.subplots(2,2, figsize=(20,20)) # total average
fig2, axes2 = plt.subplots(2,2, figsize=(20,20)) # discrete blocks
fig3, axes3 = plt.subplots(2,2, figsize=(20,20)) # block averaging

def plot_rdf(midpoints, data, label, axes):
    """
    data is a list of y values for the plotting
    """
    # Plot in 2d grid
    for i, y in enumerate(data):
        row = math.floor(i / 2)
        col = i % 2
        axes[row,col].plot(midpoints, y, label=label)



#==================================================
# Average
#==================================================
#                                             cutoff, s1, s2
cc_array_2d,   midpoints = rdf.get_rdf(atoms_list, 6,  6,  6, verbose=verbose, save=True, comment=comment)
sisi_array_2d, midpoints = rdf.get_rdf(atoms_list, 6, 14, 14, verbose=verbose, save=True, comment=comment)
sic_array_2d,  midpoints = rdf.get_rdf(atoms_list, 6, 14,  6, verbose=verbose, save=True, comment=comment)
all_array_2d,  midpoints = rdf.get_rdf(atoms_list, 6,         verbose=verbose, save=True, comment=comment)

data = [np.mean(cc_array_2d,   axis=0),
        np.mean(sisi_array_2d, axis=0),
        np.mean(sic_array_2d,  axis=0),
        np.mean(all_array_2d,  axis=0)]

plot_rdf(midpoints, data, comment, axes1)


#==================================================
# Discrete Blocks
#==================================================
n_blocks      = 20
block_indexes = np.linspace(0, 1, n_blocks+1) * all_array_2d.shape[0]
block_indexes = block_indexes.astype(int)

for k in range(n_blocks):
    data = [np.mean(cc_array_2d[block_indexes[k]:block_indexes[k+1],:],   axis=0),
            np.mean(sisi_array_2d[block_indexes[k]:block_indexes[k+1],:], axis=0),
            np.mean(sic_array_2d[block_indexes[k]:block_indexes[k+1],:],  axis=0),
            np.mean(all_array_2d[block_indexes[k]:block_indexes[k+1],:],  axis=0)]

    plot_rdf(midpoints, data, "{} (Block: {})".format(comment, k+1), axes2)

#==================================================
# Block Averaging
#==================================================
#                                                          cutoff, blocks, s1, s2
cc_array_blocks_2d,   midpoints = rdf.block_average(atoms_list, 6,     20,  6,  6, verbose=verbose, save=True, comment=comment)
sisi_array_blocks_2d, midpoints = rdf.block_average(atoms_list, 6,     20, 14, 14, verbose=verbose, save=True, comment=comment)
sic_array_blocks_2d,  midpoints = rdf.block_average(atoms_list, 6,     20, 14,  6, verbose=verbose, save=True, comment=comment)
all_array_blocks_2d,  midpoints = rdf.block_average(atoms_list, 6,     20,         verbose=verbose, save=True, comment=comment)

# Plot in 2d grid the block averaging
for k in range(cc_array_blocks_2d.shape[0]):
    data = [cc_array_blocks_2d[k,:],
            sisi_array_blocks_2d[k,:],
            sic_array_blocks_2d[k,:],
            all_array_blocks_2d[k,:]]

    plot_rdf(midpoints, data, "{} (Block: {})".format(comment, k+1), axes3)



#==================================================
# NEW DFT
#==================================================

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
    axes2[i,j].plot(midpoints, y, label=label, color="k")
    axes3[i,j].plot(midpoints, y, label=label, color="k")


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

for axes in [axes1, axes2, axes3]:
    for (i,j), title in titles:
        axes[i,j].set_title("RDF for {} bonds".format(title), fontsize=fontsize+2)

for ax in list(axes1.flatten()) + list(axes2.flatten()) + list(axes3.flatten()):
    ax.tick_params(axis="both", which="major", labelsize=fontsize-2)
    ax.set_xlabel("Distance [$\AA$]", fontsize=fontsize)
    ax.set_ylabel("g(r)",             fontsize=fontsize)
    ax.legend(ncol=3, bbox_to_anchor=(0.5, -0.1), loc="upper center", fontsize=fontsize-4)
    ax.grid()

for fig in [fig1, fig2, fig3]:
    fig.tight_layout(pad=3)

#==================================================
# Saving
#==================================================
fig1.savefig("rdf_{}.png".format(comment),                 dpi=300, bbox_inches="tight")
fig2.savefig("rdf_discrete-blocks_{}.png".format(comment), dpi=300, bbox_inches="tight")
fig3.savefig("rdf_block-averaging_{}.png".format(comment), dpi=300, bbox_inches="tight")

print("Finished!")
