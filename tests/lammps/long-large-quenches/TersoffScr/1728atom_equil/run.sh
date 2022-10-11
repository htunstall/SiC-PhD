#!/bin/bash

ulimit -s unlimited
ulimit -v unlimited



module load GCC/8.3.0 Python/3.7.4 OpenBLAS/0.3.7 OpenMPI/3.1.4 ScaLAPACK/2.0.2 netCDF/4.7.1

export OMP_NUM_THREADS=1
exe="/storage/mssnkt_grp/msufgx/gits/lammps-29Oct20_atomistica/src/lmp_mpi"




# Run the code
mpirun -np 28  ${exe} -in in.* -log lammps.log> lammps.std.out 2> lammps.std.err
