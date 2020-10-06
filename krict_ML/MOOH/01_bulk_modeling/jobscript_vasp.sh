#!/bin/sh
#PBS -V
#PBS -A vasp
#PBS -N 29_HgHO2_np
#PBS -q normal
#PBS -l select=4:ncpus=68:mpiprocs=64:ompthreads=1
#PBS -l walltime=06:00:00
#PBS -m abe
#PBS -W sandbox=PRIVATE

cd $PBS_O_WORKDIR

source /home01/x1888a07/env_vasp.sh

# rm ./POTCAR

#cat $POT/H/POTCAR >> ./POTCAR       # H   250 eV 1
#cat $POT/Ni_pv/POTCAR >> ./POTCAR   # Ni_pv 368 eV 16
#cat $POT/O/POTCAR >> ./POTCAR       # O 400 eV 6

mpirun $VASP_STD > stdout
