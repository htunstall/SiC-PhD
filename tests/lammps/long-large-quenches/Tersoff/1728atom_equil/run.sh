#!/bin/bash

ulimit -s unlimited
ulimit -v unlimited

module load GCC/6.4.0-2.28 OpenMPI/2.1.1 OpenBLAS/0.2.20 ScaLAPACK/2.0.2-OpenBLAS-0.2.20 FFTW/3.3.6 netCDF/4.4.1.1 Python/3.6.3 matplotlib/2.0.2-Python-3.6.3 ASE/3.16.2-Python-3.6.3

export OMP_NUM_THREADS=1
exe="/storage/mssnkt_grp/msufgx/gits/lammps-29Oct20/src/lmp_mpi"




# Run the code
mpirun -np 28  ${exe} -in in.* -log lammps.log> lammps.std.out 2> lammps.std.err
