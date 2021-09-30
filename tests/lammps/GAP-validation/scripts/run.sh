#!/bin/bash

ulimit -s unlimited
ulimit -v unlimited



module r QUIP_v3-2
export OMP_NUM_THREADS=1
exe="/storage/chem/msufgx/postgrad/software/lammps-29Oct20/src/lmp_mpi"



# Run the code
# srun
mpirun -np 27  ${exe} -in in.* -log lammps.log> lammps.std.out 2> lammps.std.err

#----------------
# Analysis
#----------------
# Liquid
if [ -f dump.3000K-equil.nc ]; then
    python ../../scripts/matscipy-rings.py dump.3000K-equil.nc 10 
    python ../../scripts/get_test_rdf.py dump.3000K-equil.nc 100 npt 
fi

# Quench
if [ -f dump.300K-equil.nc ]; then 
    python ../../scripts/matscipy-rings.py dump.300K-equil.nc 10 
    python ../../scripts/get_test_rdf.py dump.300K-equil.nc 100 amorphous True
fi
