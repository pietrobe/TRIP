#!/bin/bash -l

#SBATCH --ntasks=16384
#SBATCH --cpus-per-task=1
#SBATCH --constrain=gpu
#SBATCH --time=08:00:00
#SBATCH --account=u0
#SBATCH --job-name="TRIP"

srun ./solar_3D -ksp_monitor -ksp_view -ksp_rtol 1e-8



