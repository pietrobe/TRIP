#!/bin/bash -l

#SBATCH --ntasks=12288
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --account=iac90
#SBATCH --job-name="PRD_3D"
#SBATCH --qos=gp_resa
#SBATCH --nodes=110
#SBATCH --exclusive    

srun solar_3D -ksp_monitor -ksp_monitor_true_residual  -ksp_view -ksp_rtol 1e-4 -ksp_max_it 1
 
