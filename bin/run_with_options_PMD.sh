#!/bin/bash -l

### to be defined, 12288 is tested !!!!!!!!!
### #SBATCH --ntasks=12288
### #SBATCH --ntasks=16384
#SBATCH --ntasks=512

## to be defined !!!!!!!!! 
## #SBATCH --time=02:15:00

#### #SBATCH --cpus-per-task=1
#### #SBATCH --account=iac90
#SBATCH --job-name="TRIP_PRD_3D"
### #SBATCH --qos=gp_resa

#SBATCH --mail-type=ALL
#SBATCH --mail-user=simone.riva@usi.ch

#SBATCH --exclusive
##### # SBATCH --mem_per_cpu=
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --account=u2
#SBATCH --partition=debug
#SBATCH --ntasks-per-socket=46

### #SBATCH --job-name="TRIP"

export CRD=" --CRD "
export CRD="  "

export MPICH_GPU_SUPPORT_ENABLED=1
ulimit -c 0
ulimit -l unlimited
ulimit -a

# Set core dump path
# export COREDUMP_DIR=/capstor/scratch/cscs/sriva/core_dumps
# mkdir -p $COREDUMP_DIR
# echo "$COREDUMP_DIR/core.%e.%p.%t" | sudo tee /proc/sys/kernel/core_pattern
# ulimit -c unlimited

### CSCS

### export MAIN_DATA_DIR=/gpfs/projects/iac90/PORTA
export MAIN_DATA_DIR=/capstor/scratch/cscs/sriva/PORTA

# export MAIN_DATA_DIR=/gpfs/projects/iac90/input/DataSet_TRIP_PORTA/Input_Data_TRIP_PORTA_64x64
# export SERIES_DIR=Comparison-TRIP-PORTA-20240827T115438Z-001-CRD
# export GRID_SIZE=64

### MareNostrum5
#export MAIN_DATA_DIR=/gpfs/projects/iac90
#export SERIES_DIR=Comparison-TRIP-PORTA-20240827T115438Z-001-CRD
#export GRID_SIZE=32

export INPUT_DIR=${MAIN_DATA_DIR}/
export OUTPUT_DIR=${MAIN_DATA_DIR}/output_PRD_T512_B/

## For 32 x 32 grid size
export INPUT_PMD=cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd

# export INPUT_CONFIG=AR_385_Cut_32x32-CRD_I_B_V0.conf
# export INPUT_CONFIG=AR_385_Cut_32x32-CRD_I_B_V.conf
# export INPUT_CONFIG=AR_385_Cut_32x32-CRD_I_B0_V.conf

export RTOL=1e-9
export GMRES_RESTART0=30
export MAX_ITERATIONS=15

echo "INPUT_DIR: $INPUT_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "INPUT_CONFIG: $INPUT_CONFIG"
echo "CRD: $CRD"

export SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "srun --cpu-bind=socket
        ${PWD}/mps-wrapper.sh
        ${APP_PATH}/solar_3D $CRD
        --input_dir $INPUT_DIR
        --problem_pmd_file $INPUT_PMD
        --output_dir $OUTPUT_DIR
        -ksp_type fgmres
        -ksp_gmres_restart $GMRES_RESTART0
        -ksp_max_it $MAX_ITERATIONS
        -ksp_monitor
        -ksp_view
        -ksp_rtol $RTOL"

echo ""
echo "Starting TRIP ...... "
echo ""


export APP_PATH=/users/sriva/git/solar_3d/build


srun --cpu-bind=socket   \
      ${PWD}/mps-wrapper.sh \
      ${APP_PATH}/solar_3D $CRD \
      --input_dir $INPUT_DIR  \
      --problem_pmd_file $INPUT_PMD \
      --output_dir $OUTPUT_DIR \
      -ksp_type fgmres \
      -ksp_gmres_restart $GMRES_RESTART0 \
      -ksp_max_it $MAX_ITERATIONS \
      -ksp_monitor \
      -ksp_view \
      -ksp_rtol $RTOL
