#!/bin/bash -l

### to be defined, 12288 is tested !!!!!!!!!
### #SBATCH --ntasks=12288
### #SBATCH --ntasks=16384
#SBATCH --ntasks=2048

## to be defined !!!!!!!!! 
## #SBATCH --time=02:15:00

#### #SBATCH --cpus-per-task=1
#### #SBATCH --account=iac90
#SBATCH --job-name="TRIP_PRD_3D"
### #SBATCH --qos=gp_resa

#SBATCH --mail-type=ALL
#SBATCH --mail-user=simone.riva@usi.ch

#SBATCH --cpus-per-task=1
#SBATCH --constrain=gpu
#SBATCH --time=00:15:00
#SBATCH --account=u0
### #SBATCH --job-name="TRIP"

# export CRD="--CRD"
export CRD=" "

### CSCS

### export MAIN_DATA_DIR=/gpfs/projects/iac90/PORTA
export MAIN_DATA_DIR=/scratch/snx3000/sriva/PORTA

# export MAIN_DATA_DIR=/gpfs/projects/iac90/input/DataSet_TRIP_PORTA/Input_Data_TRIP_PORTA_64x64
# export SERIES_DIR=Comparison-TRIP-PORTA-20240827T115438Z-001-CRD
# export GRID_SIZE=64

### MareNostrum5
#export MAIN_DATA_DIR=/gpfs/projects/iac90
#export SERIES_DIR=Comparison-TRIP-PORTA-20240827T115438Z-001-CRD
#export GRID_SIZE=32

export INPUT_DIR=${MAIN_DATA_DIR}/
export OUTPUT_DIR=${MAIN_DATA_DIR}/output_test_extra_mu_XI_dmu/

## For 32 x 32 grid size
export INPUT_PMD=cai_0Bx_0By_0Bz_1Vx_1Vy_1Vz_GT4_5x5x133_it100.pmd

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


export APP_PATH=/users/sriva/git/solar_3d/build


srun ${APP_PATH}/solar_3D $CRD --input_dir $INPUT_DIR  --problem_pmd_file $INPUT_PMD --output_dir $OUTPUT_DIR -ksp_type fgmres -ksp_gmres_restart 30 -ksp_max_it 14 -ksp_monitor -ksp_view -ksp_rtol 1e-2




