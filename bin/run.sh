#!/bin/bash -l

#SBATCH --ntasks=2048
#SBATCH --cpus-per-task=1
#SBATCH --constrain=gpu
#SBATCH --time=04:00:00
#SBATCH --account=sm74
#SBATCH --job-name="PRT_3D"

srun ./solar_3D -ksp_monitor -ksp_view -ksp_max_it 2 
#-ksp_rtol 1e-9



