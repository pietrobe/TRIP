#!/bin/bash -l

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=mc
#SBATCH --time=05:00:00
#SBATCH --account=u2
#SBATCH --job-name="PRD_3D"

srun ./solar_3D_serial 12 8 -ksp_monitor


