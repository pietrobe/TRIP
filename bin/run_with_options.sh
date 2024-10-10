#!/bin/bash -l

#SBATCH --ntasks=2048
#SBATCH --cpus-per-task=1
#SBATCH --constrain=gpu
#SBATCH --time=00:15:00
#SBATCH --account=u0
#SBATCH --job-name="TRIP"

export CRD="--CRD"

### CSCS
export MAIN_DATA_DIR=/scratch/snx3000/sriva/test_3d

### MareNostrum 5
# export MAIN_DATA_DIR=/gpfs/projects/iac90/

export SERIES_DIR=Comparison-TRIP-PORTA-20240827T115438Z-001-CRD

### Set the horizontal grid size
export GRID_SIZE=64

### MareNostrum5
#export MAIN_DATA_DIR=/gpfs/projects/iac90
#export SERIES_DIR=Comparison-TRIP-PORTA-20240827T115438Z-001-CRD
#export GRID_SIZE=32

export INPUT_DIR=${MAIN_DATA_DIR}/input/${SERIES_DIR}/a${GRID_SIZE}
export OUTPUT_DIR=${MAIN_DATA_DIR}/output_2048_BV0_mu1/${SERIES_DIR}/a${GRID_SIZE}

## For 32 x 32 grid size
# export INPUT_CONFIG=AR_385_Cut_32x32-CRD_I_B0_V0.conf

# export INPUT_CONFIG=AR_385_Cut_64x64_mirrorxy-CRD_I_B0_V0.conf

export INPUT_CONFIG=AR_385_Cut_64x64_mirrorxy-CRD_I_B0_V0_KQ.conf


# export INPUT_CONFIG=AR_385_Cut_32x32-CRD_I_B_V0.conf
# export INPUT_CONFIG=AR_385_Cut_32x32-CRD_I_B_V.conf
# export INPUT_CONFIG=AR_385_Cut_32x32-CRD_I_B0_V.conf

echo "INPUT_DIR: $INPUT_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "INPUT_CONFIG: $INPUT_CONFIG"
echo "CRD: $CRD"
echo "srun ./solar_3D $CRD --input_dir $INPUT_DIR  --problem_input_config $INPUT_CONFIG --output_dir $OUTPUT_DIR -ksp_type fgmres -ksp_gmres_restart 30 -ksp_max_it 200 -ksp_monitor -ksp_view -ksp_rtol 1e-8"

echo ""
echo "Running the simulation..."
echo ""

srun ./solar_3D $CRD --input_dir $INPUT_DIR  --problem_input_config $INPUT_CONFIG --output_dir $OUTPUT_DIR -ksp_type fgmres -ksp_gmres_restart 30 -ksp_max_it 200 -ksp_monitor -ksp_view -ksp_rtol 1e-8




