#!/bin/bash -l

#SBATCH --ntasks=12288
#SBATCH --cpus-per-task=1
#SBATCH --constrain=gpu
#SBATCH --time=00:50:00
#SBATCH --account=sm74
#SBATCH --job-name="TRIP"

srun ./solar_3D -ksp_monitor -ksp_view -ksp_rtol 1e-9



