#!/bin/bash -l

#SBATCH --ntasks=1152
#SBATCH --cpus-per-task=1
#SBATCH --constrain=gpu
#SBATCH --time=01:00:00
#SBATCH --account=sm74
#SBATCH --job-name="PRD_3D"

srun ./solar_3D -ksp_monitor -ksp_rtol 1e-9 -ksp_view



