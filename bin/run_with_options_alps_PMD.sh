#!/bin/bash -l

### to be defined, 12288 is tested !!!!!!!!!
### #SBATCH --ntasks=12288
### #SBATCH --ntasks=16384
#SBATCH --ntasks=6144

## to be defined !!!!!!!!! 
#SBATCH --time=01:30:00

#SBATCH --job-name="TRIP_PRD_3D"

## #SBATCH --mail-type=ALL
## #SBATCH --mail-user=my.email

#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --account=c40
#SBATCH --partition=normal
#SBATCH --ntasks-per-socket=48



echo ""
echo "=========================================="
echo "SLURM Environment Variables:"
echo "=========================================="
printf "%-30s %s\n" "SLURM_JOB_ID:" "$SLURM_JOB_ID"
printf "%-30s %s\n" "SLURM_JOB_NAME:" "$SLURM_JOB_NAME"
printf "%-30s %s\n" "SLURM_JOB_NODELIST:" "$SLURM_JOB_NODELIST"
printf "%-30s %s\n" "SLURM_JOB_NUM_NODES:" "$SLURM_JOB_NUM_NODES"
printf "%-30s %s\n" "SLURM_NTASKS:" "$SLURM_NTASKS"
printf "%-30s %s\n" "SLURM_NTASKS_PER_NODE:" "$SLURM_NTASKS_PER_NODE"
printf "%-30s %s\n" "SLURM_CPUS_PER_TASK:" "$SLURM_CPUS_PER_TASK"
printf "%-30s %s\n" "SLURM_CPUS_ON_NODE:" "$SLURM_CPUS_ON_NODE"
printf "%-30s %s\n" "SLURM_MEM_PER_NODE:" "$SLURM_MEM_PER_NODE"
printf "%-30s %s\n" "SLURM_SUBMIT_DIR:" "$SLURM_SUBMIT_DIR"
printf "%-30s %s\n" "SLURM_SUBMIT_HOST:" "$SLURM_SUBMIT_HOST"
printf "%-30s %s\n" "SLURM_NODEID:" "$SLURM_NODEID"
printf "%-30s %s\n" "SLURM_PROCID:" "$SLURM_PROCID"
printf "%-30s %s\n" "SLURM_LOCALID:" "$SLURM_LOCALID"
printf "%-30s %s\n" "SLURM_GPUS:" "$SLURM_GPUS"
printf "%-30s %s\n" "SLURM_GPUS_PER_NODE:" "$SLURM_GPUS_PER_NODE"
printf "%-30s %s\n" "SLURM_OUTPUT_LOG:" "slurm-${SLURM_JOB_ID}.out"
echo "=========================================="
echo ""


export CRD=" --CRD "
export CRD="  "

export MPICH_GPU_SUPPORT_ENABLED=1

ulimit -c 0
ulimit -l unlimited
ulimit -a

export MAIN_DATA_DIR=/capstor/scratch/cscs/sriva/PORTA

export INPUT_DIR=${MAIN_DATA_DIR}/
export OUTPUT_DIR=${MAIN_DATA_DIR}/output_PRD_T32x32_AR/

## For 32 x 32 grid size
export INPUT_PMD=cai_1Bx_1By_1Bz_1Vx_1Vy_1Vz_GT4_32x32x133.pmd

export RTOL=1e-5
export GMRES_RESTART=30
export MAX_ITERATIONS=2

echo "INPUT_DIR: $INPUT_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "INPUT_CONFIG: $INPUT_CONFIG"
echo "CRD: $CRD"

export APP_PATH=${HOME}/git/solar_3d/build
export SCRIPT_DIR=${HOME}/git/solar_3d/bin

echo ""
echo "=========================================="
echo "User-Defined Variables:"
echo "=========================================="
printf "%-30s %s\n" "CRD:" "$CRD"
printf "%-30s %s\n" "MPICH_GPU_SUPPORT_ENABLED:" "$MPICH_GPU_SUPPORT_ENABLED"
printf "%-30s %s\n" "MAIN_DATA_DIR:" "$MAIN_DATA_DIR"
printf "%-30s %s\n" "INPUT_DIR:" "$INPUT_DIR"
printf "%-30s %s\n" "OUTPUT_DIR:" "$OUTPUT_DIR"
printf "%-30s %s\n" "INPUT_PMD:" "$INPUT_PMD"
printf "%-30s %s\n" "RTOL:" "$RTOL"
printf "%-30s %s\n" "GMRES_RESTART:" "$GMRES_RESTART"
printf "%-30s %s\n" "MAX_ITERATIONS:" "$MAX_ITERATIONS"
printf "%-30s %s\n" "APP_PATH:" "$APP_PATH"
printf "%-30s %s\n" "SCRIPT_DIR:" "$SCRIPT_DIR"
echo "=========================================="
echo ""


# Define arguments as a bash array
ARGS=(
    "$CRD"
    "--input_dir" "$INPUT_DIR"
    "--problem_pmd_file" "$INPUT_PMD"
    "--output_dir" "$OUTPUT_DIR"
    "-ksp_type" "fgmres"
    "-ksp_gmres_restart" "$GMRES_RESTART"
    "-ksp_max_it" "$MAX_ITERATIONS"
    "-ksp_monitor"
    "-ksp_view"
    "-ksp_rtol" "$RTOL"
)

echo "Running command with MPS wrapper:"
echo " srun --cpu-bind=socket ${SCRIPT_DIR}/mps-wrapper.sh ${APP_PATH}/solar_3D ${ARGS[@]}"
echo ""
echo ""
echo "Starting TRIP ...... "
echo ""

srun --cpu-bind=socket \
      ${SCRIPT_DIR}/mps-wrapper.sh \
      ${APP_PATH}/solar_3D "${ARGS[@]}"
