#!/bin/bash -l

#SBATCH --ntasks=512
#SBATCH --nodes=64
#SBATCH --cpus-per-task=1
#SBATCH --constraint=mc
#SBATCH --time=01:10:00
#SBATCH --account=u2

#export LD_LIBRARY_PATH=/users/pietrob/rii/rii-c/build/
export LD_LIBRARY_PATH=$TRILINOS_DIR/lib:$LD_LIBRARY_PATH

srun ./main  



