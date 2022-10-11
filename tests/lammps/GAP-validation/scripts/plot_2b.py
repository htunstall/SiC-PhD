import os
import sys
import quippy
import ase.io
from   tqdm              import tqdm
import numpy             as     np
import matplotlib.pyplot as     plt
#quippy.system_module.verbosity_push(10) #quippy.PRINT_NERD)

comment = os.path.split(os.path.split(os.getcwd())[0])[1]

bonds         = ["SiSi", "SiC", "CC"]
data_per_bond = {bond: [] for bond in bonds}

if not os.path.isfile("dimer_SiSi.txt"):
    print("No saved dimer data, rerunning")
    sys.stdout.write("="*60+"\n"+"="*60+"\n")
    print("gp.xml")
    sys.stdout.flush()
    calc           = quippy.potential.Potential("IP GAP", param_filename="gap/gp.xml")
#    calc_2b_hack   = calc
#    calc_soap_hack = calc
    sys.stdout.write("="*60+"\n"+"="*60+"\n")
    print("gp-2b.xml")
    sys.stdout.flush()
    calc_2b_hack   = quippy.potential.Potential("IP GAP", param_filename="gap/gp-2b.xml")
    sys.stdout.write("="*60+"\n"+"="*60+"\n")
    print("gp-soap.xml")
    sys.stdout.flush()
    calc_soap_hack = quippy.potential.Potential("IP GAP", param_filename="gap/gp-soap.xml")
    sys.stdout.write("="*60+"\n"+"="*60+"\n")
    sys.stdout.flush()


    atoms_list     = ase.io.read("../../scripts/2b-plus.xyz", index=":")

    for a in tqdm(atoms_list):
        bond = None
        if a[0].symbol != a[1].symbol:
            bond = "SiC"
        elif a[0].symbol == "Si":
            bond = "SiSi"
        elif a[0].symbol == "C":
            bond = "CC"
        else:
            sys.exit("Can't work out the bond")
   
        dft_energy     = np.nan
        if "dft_energy" in a.info.keys():
            dft_energy = a.info["dft_energy"]

        data_per_bond[bond].append([a.info["bond_length"],                  # Distance
                                   dft_energy,                              # DFT Energy
                                   calc.get_potential_energy(a),            # GAP Energy ALL
                                   calc_2b_hack.get_potential_energy(a),    # GAP Energy 2b ONLY (HACK)
                                   calc_soap_hack.get_potential_energy(a)]) # GAP Energy SOAP ONLY (HACK)
    
    for bond in bonds:
        data_per_bond[bond] = np.array(data_per_bond[bond], dtype=float)
#        data_per_bond[bond][:,1:] -= np.nanmin(data_per_bond[bond][:,1:])
        np.savetxt("dimer_{}.txt".format(bond), data_per_bond[bond])
        
else:
    print("Loading dimer data from file")
    for bond in bonds:
        data_per_bond[bond] = np.loadtxt("dimer_{}.txt".format(bond))


print("Plotting")
fontsize=20
fig, axes  = plt.subplots(1, 3, figsize=(30,10))

for ax, bond in zip(axes, bonds):
    ax.set_title("Bond {}".format(bond), fontsize=fontsize+4)
    data = data_per_bond[bond]
    ax.plot(data[:,0], data[:,1]/2,  "-o", c="k", label="DFT")
    ax.plot(data[:,0]+0.01, data[:,2]/2,  "-o",        label="{} (ALL)".format(comment))
    ax.plot(data[:,0]+0.02, data[:,3]/2,  "-o",        label="{} (2b Hack)".format(comment))
    ax.plot(data[:,0]+0.03, data[:,4]/2,  "-o",        label="{} (SOAP Hack)".format(comment))

    ax.grid()
    ax.tick_params(axis="both", labelsize=fontsize)


for _ax in axes:
    _ax.set_xlabel("Distance [$\AA$]", fontsize=fontsize)
    _ax.set_ylabel("Energy per Atom [eV/atom]", fontsize=fontsize)
    top = min(50, _ax.get_ylim()[1])
    _ax.set_ylim(top=top)

axes[0].set_xlim(0, 4.5)
axes[1].set_xlim(0, 4.0)
axes[2].set_xlim(0, 4.0)


handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.025), ncol=1, fontsize=fontsize)

fig.savefig("dimer-energies.png", dpi=200, bbox_inches="tight")
