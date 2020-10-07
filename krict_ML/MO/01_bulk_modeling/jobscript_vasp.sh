#!/bin/sh
#PBS -V
#PBS -N 29_HgO_np
#PBS -q normal
#PBS -A vasp
#PBS -l select=4:ncpus=68:mpiprocs=64:ompthreads=1
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -W sandbox=PRIVATE

cd $PBS_O_WORKDIR

module load craype-mic-knl intel/18.0.3 impi/18.0.3
source /home01/x1895a04/vasp/vasp.5.4.4/env_vasp.sh

#rm ./POTCAR
# cat $POT/O/POTCAR > ./POTCAR      # O    520 eV 4
# cat $POT/V_pv/POTCAR >> ./POTCAR     # V    520 eV 4
# cat $POT/O/POTCAR >> ./POTCAR      # O     400 eV 6  
 
mpirun $VASP/vasp_std > stdout

