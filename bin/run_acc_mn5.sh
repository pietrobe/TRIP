#!/bin/bash -l

#### #SBATCH --ntasks=6912
#SBATCH --ntasks=512

## to be defined !!!!!!!!! 
#SBATCH --time=00:59:59

#SBATCH --cpus-per-task=1
#SBATCH --account=ehpc238
#SBATCH --job-name="TRIP_PRD_3D"
#SBATCH --qos=acc_debug

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

export APP_PATH=${HOME}/git/TRIP/build
export SCRIPT_DIR=${HOME}/git/TRIP/bin

echo ""
echo "=========================================="
echo "User-Defined Variables:"
echo "=========================================="
printf "%-30s %s\n" "APP_PATH:" "$APP_PATH"
printf "%-30s %s\n" "SCRIPT_DIR:" "$SCRIPT_DIR"
echo "=========================================="
echo ""

echo ""
echo "Starting TRIP ...... "
echo ""

export OMPI_MCA_mpi_cuda_support=1

# Define arguments as a bash array
ARGS=(
    "--yaml_config" "$HOME/git/TRIP/bin/config_acc_mn5.yml" 
    "-ksp_monitor"
    "-ksp_view"
)

echo "Executing command:"
echo " srun --cpu-bind=socket  mpirun --bind-to none  ${SCRIPT_DIR}/mnacc_wrapper.sh ${APP_PATH}/solar_3D ${ARGS[@]}"
echo ""

mpirun --bind-to none  \
      ${SCRIPT_DIR}/mnacc_wrapper.sh \
      ${APP_PATH}/solar_3D "${ARGS[@]}"
