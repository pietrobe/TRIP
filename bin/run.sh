#!/bin/bash -l

#SBATCH --ntasks=140
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --constraint=mc
#SBATCH --time=01:10:00
#SBATCH --account=u2
#SBATCH --job-name="PRD_3D"

srun ./solar_3D 4 4


