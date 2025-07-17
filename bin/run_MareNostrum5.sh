#!/bin/bash -l

#SBATCH --ntasks=512

#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --account=ehpc238
#SBATCH --job-name="TRIP_PRD_3D"
#SBATCH --qos=gp_debug
#SBATCH --exclusive

srun ./solar_3D -ksp_monitor -ksp_view -ksp_rtol 1e-5 #-ksp_max_it 1 #-info -help -log_view
