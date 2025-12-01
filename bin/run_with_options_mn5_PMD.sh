#!/bin/bash -l

#SBATCH --ntasks=512

## to be defined !!!!!!!!! 
#SBATCH --time=00:36:27

#SBATCH --cpus-per-task=1
#SBATCH --account=ehpc238
#SBATCH --job-name="TRIP_PRD_3D"
#SBATCH --qos=acc_debug

#SBATCH --mail-type=ALL
#SBATCH --mail-user=my.email

#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4

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

export OMPI_MCA_opal_cuda_support=1
export MAIN_DATA_DIR=/gpfs/projects/ehpc238/PORTA

export INPUT_DIR=${MAIN_DATA_DIR}/
export OUTPUT_DIR=${MAIN_DATA_DIR}/output_PRD_T512_B/

## For 32 x 32 grid size
export INPUT_PMD=cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd

export RTOL=1e-6
export GMRES_RESTART=30
export MAX_ITERATIONS=15

export APP_PATH=${HOME}/git/solar_3d/build
export SCRIPT_DIR=${HOME}/git/solar_3d/bin

echo ""
echo "=========================================="
echo "User-Defined Variables:"
echo "=========================================="
printf "%-30s %s\n" "CRD:" "$CRD"
printf "%-30s %s\n" "OMPI_MCA_opal_cuda_support:" "$OMPI_MCA_opal_cuda_support"
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

echo ""
echo "Starting TRIP ...... "
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

echo "Executing command:"
echo " srun --cpu-bind=socket  mpirun --bind-to none  ${SCRIPT_DIR}/mnacc_wrapper.sh ${APP_PATH}/solar_3D ${ARGS[@]}"
echo ""


mpirun --bind-to none  \
      ${SCRIPT_DIR}/mnacc_wrapper.sh \
      ${APP_PATH}/solar_3D "${ARGS[@]}"
