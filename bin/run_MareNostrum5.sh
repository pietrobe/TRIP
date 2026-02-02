#!/bin/bash -l

#SBATCH --ntasks=6144
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --account=ehpc238
#SBATCH --job-name="TRIP_PRD_3D"
#SBATCH --qos=gp_ehpc
#SBATCH --exclusive

srun ./solar_3D
