#!/bin/bash -l


#SBATCH --ntasks=2048
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=1
#SBATCH --time=01:00:00
#SBATCH --account=u0
#SBATCH --job-name="TRIP"
#SBATCH --constrain=gpu
#SBATCH --partition=normal
#

export OUT_DATA_DIR=/scratch/snx3000/sriva/test_3d
export IN_DATA_DIR=${HOME}/git/solar_3d/input/PORTA
export ATHMOSPHERE_NAME=cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd
	
srun  ${HOME}/git/solar_3d/build/solar_3D --CRD --input_dir ${IN_DATA_DIR} --output_dir ${OUT_DATA_DIR} --problem_input_file  ${ATHMOSPHERE_NAME} -ksp_monitor -ksp_monitor_true_residual  -ksp_view -ksp_rtol 1e-8 -ksp_max_it 130


