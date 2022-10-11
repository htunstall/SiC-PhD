import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams['agg.path.chunksize'] = 10000

filename = "lammps.std.out"

with open(filename, "r") as f:
    lines = f.readlines()

start_is     = np.array([i for i, line in enumerate(lines) if line.startswith("Per MPI rank memory allocation")]) + 2
data_lengths = []
j = 0
for i, line in enumerate(lines):
    if line.startswith("Loop time of "):
        data_lengths.append(int(i - start_is[j]))
        j += 1

if len(data_lengths) == 0:
    data_lengths.append(int(i - start_is[j] - 1))

data = np.loadtxt(filename, skiprows=start_is[0], max_rows=data_lengths[0])
for j, i in enumerate(start_is[1:]):
    data = np.concatenate((data, np.loadtxt(filename, skiprows=i, max_rows=data_lengths[j+1])), axis=0)


# 0    1    2    3      4      5      6     7       8      9  10 11 12  13  14  15  16  17
# Step Time Temp PotEng KinEng TotEng Press Density Volume Lx Ly Lz Xlo Xhi Ylo Yhi Zlo Zhi

data[:,6] = data[:,6] * 0.0001 # bar -> GPa

step = data[:,0]
time = data[:,1] # ps
temp = data[:,2] # K
pote = data[:,3] # eV
kine = data[:,4] # eV
tote = data[:,5] # eV
pres = data[:,6] # GPa
dens = data[:,7] # g/cm^3
vol  = data[:,8] # Angs^3

fig, axes = plt.subplots(2,4, figsize=(40,20))

# kin
axes[0,0].plot(time, kine - np.amin(kine))
axes[0,0].set_ylabel("Energy [eV]", fontsize=16)
axes[0,0].set_title("Kinetic Energy", fontsize=18)

# pot
axes[0,1].plot(time, pote - np.amin(pote))
axes[0,1].set_ylabel("Energy [eV]", fontsize=16)
axes[0,1].set_title("Potential Energy", fontsize=18)

# tot
axes[0,2].plot(time, tote - np.amin(tote))
axes[0,2].set_ylabel("Energy [eV]", fontsize=16)
axes[0,2].set_title("Total Energy", fontsize=18)

# temp
axes[1,0].plot(time, temp)
axes[1,0].set_ylabel("Temperature [K]", fontsize=16)
axes[1,0].set_title("Temperature", fontsize=18)

# pres
axes[1,1].plot(time, pres)
axes[1,1].set_ylabel("Pressure [GPa]", fontsize=16)
axes[1,1].set_title("Pressure", fontsize=18)

# vol
axes[1,2].plot(time, vol)
axes[1,2].set_ylabel("Volume [$\AA^3$]", fontsize=16)
axes[1,2].set_title("Volume", fontsize=18)

# dens
axes[1,3].plot(time, dens)
axes[1,3].set_ylabel("Desnity [$g/cm^3$]", fontsize=16)
axes[1,3].set_title("Density", fontsize=18)

for ax in axes.flatten():
    ax.set_xlabel("Time [ps]", fontsize=16)
    ax.tick_params(axis="both", which="both", labelsize=16)
    ax.grid()
#    ax.set_xlim(0,10)

fig.tight_layout()
fig.savefig("physical-quantaties.png", dpi=300)

np.savetxt("physical-data.txt", data)
