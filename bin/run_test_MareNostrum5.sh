#!/bin/bash -l

## to be defined !!!!!!!!!
#SBATCH --ntasks=1024

#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --account=iac90
#SBATCH --job-name="TRIP_PRD_3D"
#SBATCH --qos=gp_resa
## #SBATCH --nodes=1

### set here the main input data directory
export MAIN_INPUT_PROBLEMS_SERIES_DIR=PORTA

### set here the input and output data directories
export IN_DATA_DIR=/home/usi/usi441290/git/solar_3d/input/${MAIN_INPUT_PROBLEMS_SERIES_DIR}

### the main output directory
export OUT_DATA_DIR=/gpfs/projects/iac90/output_test/${MAIN_INPUT_PROBLEMS_SERIES_DIR}

### set here the atmosphere name
export ATHMOSPHERE_NAME=cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd

### ATTENTION: remove --CRD to perform a PRD calculation
srun ./solar_3D --CRD --input_dir ${IN_DATA_DIR} --output_dir ${OUT_DATA_DIR} --problem_input_file  ${ATHMOSPHERE_NAME} -ksp_monitor -ksp_monitor_true_residual  -ksp_view -ksp_rtol 1e-6  -ksp_max_it 130
