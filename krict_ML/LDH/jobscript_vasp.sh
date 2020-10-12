#!/bin/sh
#PBS -V
#PBS -N MnO2H2
#PBS -q normal
#PBS -l select=4:ncpus=64:mpiprocs=64
#PBS -l walltime=06:00:00
#PBS -m abe
#PBS -A vasp
#PBS -W sandbox=PRIVATE

cd $PBS_O_WORKDIR

source /home01/x2045a06/vasp/env_vasp.sh
 
mpirun $VASP_STD > stdout

