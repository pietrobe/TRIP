#!/bin/bash -l

#SBATCH --ntasks=6912

#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --account=iac90
#SBATCH --job-name="TRIP_PRD_3D"
#SBATCH --qos=gp_resa
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96

srun ./solar_3D -ksp_monitor -ksp_view -ksp_rtol 1e-8
