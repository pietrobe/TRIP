#!/bin/bash -l

#SBATCH --ntasks=2048
#SBATCH --cpus-per-task=1
#SBATCH --constrain=mc
#SBATCH --time=02:00:00
#SBATCH --account=u0
#SBATCH --job-name="TRIP"

export CRD="--CRD"

### CSCS
export MAIN_DATA_DIR=/scratch/snx3000/sriva/test_3d/input/orig
export SERIES_DIR=Comparison-TRIP-PORTA-20240827T115438Z-001-CRD
export GRID_SIZE=32

### MareNostrum5
#export MAIN_DATA_DIR=/gpfs/projects/iac90
#export SERIES_DIR=Comparison-TRIP-PORTA-20240827T115438Z-001-CRD
#export GRID_SIZE=32

export INPUT_DIR=${MAIN_DATA_DIR}/
export OUTPUT_DIR=${MAIN_DATA_DIR}/output/

## For 32 x 32 grid size
export INPUT_PMD=cai_1Bx_1By_1Bz_1Vx_1Vy_1Vz_GT4_32x32x133.pmd
# export INPUT_CONFIG=AR_385_Cut_32x32-CRD_I_B_V0.conf
# export INPUT_CONFIG=AR_385_Cut_32x32-CRD_I_B_V.conf
# export INPUT_CONFIG=AR_385_Cut_32x32-CRD_I_B0_V.conf

echo "INPUT_DIR: $INPUT_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "INPUT_CONFIG: $INPUT_CONFIG"
echo "CRD: $CRD"
echo "srun ./solar_3D $CRD --input_dir $INPUT_DIR  --problem_pmd_file $INPUT_PMD --output_dir $OUTPUT_DIR -ksp_type fgmres -ksp_gmres_restart 30 -ksp_max_it 200 -ksp_monitor -ksp_view -ksp_rtol 1e-8"

echo ""
echo "Running the simulation..."
echo ""

srun ./solar_3D $CRD --input_dir $INPUT_DIR  --problem_pmd_file $INPUT_PMD --output_dir $OUTPUT_DIR -ksp_type fgmres -ksp_gmres_restart 30 -ksp_max_it 200 -ksp_monitor -ksp_view -ksp_rtol 1e-8




