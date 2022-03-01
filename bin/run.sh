#!/bin/bash -l

#SBATCH --ntasks=70
#SBATCH --nodes=5
#SBATCH --cpus-per-task=1
#SBATCH --constraint=mc
#SBATCH --time=01:10:00
#SBATCH --account=u2
#SBATCH --job-name="PRD_3D"

srun ./solar_3D 2 4 -ksp_monitor


