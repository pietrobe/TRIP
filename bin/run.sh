#!/bin/bash -l

#SBATCH --ntasks=70
#SBATCH --nodes=4
#SBATCH --cpus-per-task=1
#SBATCH --constraint=mc
#SBATCH --time=03:00:00
#SBATCH --account=u2
#SBATCH --job-name="PRD_3D"

srun ./solar_3D_test 12 8 -ksp_monitor



