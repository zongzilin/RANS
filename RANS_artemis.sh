#!/bin/bash
#PBS -P COMPENG
#PBS -q defaultQ 
#PBS -l select=1:ncpus=24:mpiprocs=24:mem=120Gb

#PBS -l walltime=1:00:00

cd $PBS_O_WORKDIR

module load gcc
module load openmpi-gcc/4.1.1

mpirun -n 24 /project/COMPENG/zzon4574 a.out > output

 exit 0

