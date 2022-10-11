#!/bin/bash

module load GCC/6.4.0-2.28 OpenMPI/2.1.1 OpenBLAS/0.2.20 ScaLAPACK/2.0.2-OpenBLAS-0.2.20 FFTW/3.3.6 netCDF/4.4.1.1 Python/3.6.3 matplotlib/2.0.2-Python-3.6.3 ASE/3.16.2-Python-3.6.3

exe="/storage/mssnkt_grp/msufgx/gits/lammps-29Oct20/src/lmp_mpi"

mpirun -n 4 ${exe} -in in.* > lammps.std.out 2> lammps.std.err

python plot-time-mem.py

python plot-lammps-physical-quantities.py

mkdir 300K
cd 300K
ln -s ../dump.300K-equil.nc dump.nc

cd ..
