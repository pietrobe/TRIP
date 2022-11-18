#!/bin/bash -l

#SBATCH --ntasks=96
#SBATCH --cpus-per-task=1
#SBATCH --constraint=mc
#SBATCH --time=00:30:00
#SBATCH --account=u2
#SBATCH --job-name="PRD_3D"

srun ./solar_3D 12 8 -ksp_monitor -ksp_view



