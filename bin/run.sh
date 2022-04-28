#!/bin/bash -l

#SBATCH --ntasks=288
#SBATCH --nodes=20
#SBATCH --cpus-per-task=1
#SBATCH --constraint=mc
#SBATCH --time=01:30:00
#SBATCH --account=u2
#SBATCH --job-name="PRD_3D"

export OMP_NUM_THREADS=1

srun ./solar_3D 12 8 -ksp_monitor



