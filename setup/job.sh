#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=75
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH -p cpuonly
#SBATCH -J 02
#SBATCH --output="out_%j.sh"

module purge
module load compiler/gnu/10
module load mpi/openmpi/4.1
module load lib/hdf5/1.12



. ap.setenv
make cleanrun
