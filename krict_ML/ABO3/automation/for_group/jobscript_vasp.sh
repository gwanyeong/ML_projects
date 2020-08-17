#!/bin/sh
#PBS -V
#PBS -A vasp
#PBS -N 
#PBS -q normal
#PBS -l select=4:ncpus=68:mpiprocs=68:ompthreads=1
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -M gyjung@unist.ac.kr
#PBS -W sandbox=PRIVATE

cd $PBS_O_WORKDIR

source /home01/x1907a08/jobscript/env_vasp3.sh

mpirun $VASP > stdout

