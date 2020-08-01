#!/bin/sh
#PBS -V
#PBS -A vasp
#PBS -N 
#PBS -q flat
#PBS -l select=8:ncpus=34:mpiprocs=34:ompthreads=1
#PBS -l walltime=04:00:00
#PBS -m abe
#PBS -M gyjung@unist.ac.kr
#PBS -W sandbox=PRIVATE

cd $PBS_O_WORKDIR

source /home01/x1907a08/jobscript/env_vasp3.sh

rm ./POTCAR

 cat $POT/C/POTCAR >> ./POTCAR        # C       400 eV 4
 cat $POT/H/POTCAR >> ./POTCAR        # H       250 eV 1
 cat $POT/N/POTCAR >> ./POTCAR        # N       400 eV 5
 cat $POT/O/POTCAR >> ./POTCAR        # O       400 eV 6
 cat $POT/S/POTCAR >> ./POTCAR        # S       259 eV 6

mpirun $VASP > stdout

