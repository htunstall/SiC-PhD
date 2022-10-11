#!/bin/bash

export OMP_NUM_THREADS=27

module restore QUIP_v3-2

exe="/storage/mssnkt_grp/msufgx/gits/lammps-29Oct20/src/lmp_mpi"

mpirun -n 4 ${exe} -in in.* > lammps.std.out 2> lammps.std.err

#python plot-time-mem.py

python plot-lammps-physical-quantities.py

mkdir 300K
cd 300K
ln -s ../dump.300K-equil.nc dump.nc

cd ..
