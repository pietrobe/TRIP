#!/bin/bash -l

## to be defined !!!!!!!!!
#SBATCH --ntasks=12288

#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --account=iac90
#SBATCH --job-name="TRIP_PRD_3D"
#SBATCH --qos=gp_resa
## #SBATCH --nodes=1

srun ./solar_3D -ksp_monitor -ksp_view -ksp_rtol 1e-6 -ksp_max_it 2
