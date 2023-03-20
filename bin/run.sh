#!/bin/bash -l

#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1
#SBATCH --constraint=mc
#SBATCH --time=00:10:00
#SBATCH --account=u2
#SBATCH --job-name="PRD_3D"

srun ./solar_3D 12 8 -ksp_monitor -ksp_view



