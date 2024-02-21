#!/bin/bash -l

#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --account=u2
#SBATCH --job-name="PRT_3D"
#SBATCH --constrain=mc

srun ./solar_3D -ksp_monitor -ksp_view -ksp_rtol 1e-9



