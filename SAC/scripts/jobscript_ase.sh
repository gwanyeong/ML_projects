#!/bin/sh
#PBS -V
#PBS -A vasp
#PBS -N dvb4_Mn
#PBS -q flat
#PBS -l select=4:ncpus=68:mpiprocs=64:ompthreads=1
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -M gyjung.jobs@gmail.com
#PBS -W sandbox=PRIVATE

cd $PBS_O_WORKDIR

/scratch/x2045a01/anaconda3/envs/cgcnn/bin/python jobs/calc_Mn.py >> logs/stdout_Mn
