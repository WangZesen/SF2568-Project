#!/bin/bash -l

# job name
#SBATCH -J myjob
# account
#SBATCH -A edu18.SF2568
# email notification
#SBATCH --mail-type=BEGIN,END
# 10 minutes wall-clock time will be given to this job
#SBATCH -t 00:02:00
# Number of nodes
#SBATCH --nodes=1
# set tasks per node to 24 in order to disablr hyperthreading
#SBATCH --ntasks-per-node=24

module add i-compilers intelmpi
mpirun -np 1 ./convex3d 10000000 0
