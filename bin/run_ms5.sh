#!/bin/bash -l

#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --account=iac90
#SBATCH --job-name="PRT_3D"
#SBATCH --qos=gp_resa
#SBATCH --nodes=1
## #SBATCH --constrain=mc

srun ./solar_3D -ksp_monitor -ksp_monitor_true_residual  -ksp_view -ksp_rtol 1e-9  -ksp_max_it 20
