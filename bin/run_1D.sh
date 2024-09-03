#!/bin/bash -l

#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=01:59:01
#SBATCH --job-name="My_PRD"

#SBATCH --partition=normal
#SBATCH --exclusive
#SBATCH --hint=multithread
#SBATCH --account=u2
#SBATCH --constraint=mc

### #SBATCH --mem=120GB

#SBATCH --mail-type=ALL
#SBATCH --mail-user=simone.rva@gmail.com

#SBATCH --output=/users/sriva/cscs_log_jobs/solar_3D_%j_log.out
#SBATCH --error=/users/sriva/cscs_log_jobs/solar_3D_%j_log.err

# export OMP_PROC_BIND=true


srun ./solar_3D  -ksp_monitor_true_residual -ksp_monitor -ksp_rtol 1e-9 -ksp_view -ksp_max_it 20 



