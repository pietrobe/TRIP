#!/bin/bash -l

#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --constrain=gpu
#SBATCH --time=00:15:00
#SBATCH --account=sm74
#SBATCH --job-name="PRD_3D"

srun ./solar_3D -ksp_monitor -ksp_view



