#!/bin/bash -l

#SBATCH --ntasks=12288
#SBATCH --nodes=96
#SBATCH --cpus-per-task=1
#SBATCH --time=01:30:00
#SBATCH --account=u2
#SBATCH --job-name="TRIP"
#SBATCH --partition=normal
#SBATCH --constrain=mc
#SBATCH --exclusive

srun solar_3D -ksp_monitor -ksp_view -ksp_rtol 1e-4  -ksp_max_it 3
 
