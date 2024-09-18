#!/bin/bash -l

## to be defined !!!!!!!!!
#SBATCH --ntasks=0

#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --account=iac90
#SBATCH --job-name="TRIP_PRD_3D"
#SBATCH --qos=gp_resa
## #SBATCH --nodes=1

## TODO in order of priority and execution.
## AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_V0_conv.pmd
## AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_conv.pmd
## AR_385_Cut_64x64_mirrorxy-CRD_I_V0-V0_conv.pmd
## AR_385_Cut_64x64_mirrorxy-CRD_I_V0_conv.pmd

### set here the main input data directory
export MAIN_INPUT_PROBLEMS_SERIES_DIR=Comparison-TRIP-PORTA-20240827T115438Z-001-CRD

### set here the input and output data directories
export IN_DATA_DIR=/gpfs/projects/iac90/input/${MAIN_INPUT_PROBLEMS_SERIES_DIR}

### the main output directory
export OUT_DATA_DIR=/gpfs/projects/iac90/output/${MAIN_INPUT_PROBLEMS_SERIES_DIR}

### set here the atmosphere name
export ATHMOSPHERE_NAME=to_be_defined.pmd

### ATTENTION: remove --CRD to perform a PRD calculation
srun ./solar_3D --CRD --input_dir ${IN_DATA_DIR} --output_dir ${OUT_DATA_DIR} --problem_input_file  ${ATHMOSPHERE_NAME} -ksp_monitor -ksp_monitor_true_residual  -ksp_view -ksp_rtol 1e-6  -ksp_max_it 130
