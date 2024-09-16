#!/bin/bash -l

#SBATCH --ntasks=2048
#SBATCH --nodes=16
#SBATCH --cpus-per-task=1
#SBATCH --time=01:30:00
#SBATCH --account=u2
#SBATCH --job-name="TRIP"
#SBATCH --constrain=mc
#SBATCH --partition=normal

srun ./solar_3D -ksp_monitor -ksp_view -ksp_rtol 1e-4



