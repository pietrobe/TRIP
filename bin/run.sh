#!/bin/bash -l

#SBATCH --ntasks=12288
#SBATCH --cpus-per-task=1
#SBATCH --constrain=gpu
#SBATCH --time=00:30:00
#SBATCH --account=u0
#SBATCH --job-name="TRIP"

srun ./solar_3D -ksp_monitor -ksp_view -ksp_rtol 1e-7



