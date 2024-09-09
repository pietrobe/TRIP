#!/bin/bash -l


#SBATCH --ntasks=2048
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=1
#SBATCH --time=04:30:00
#SBATCH --account=u0
#SBATCH --job-name="TRIP"
#SBATCH --constrain=mc
#SBATCH --partition=normal
#

srun ./solar_3D -ksp_monitor -ksp_monitor_true_residual  -ksp_view -ksp_rtol 1e-7 -ksp_max_it 20



