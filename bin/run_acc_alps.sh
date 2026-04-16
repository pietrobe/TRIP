#!/bin/bash -l

###  next: 6144 -> w_64x64x48.pmd , 12288 -> w_64x64x96.pmd  16384 -> w_64x64x128.pmd
###  next: 
###  2048 -> w_32x32x64.pmd * B
###  4096 -> w_32x32x128.pmd * B 
###  6144 -> w_64x64x48.pmd  * B 
###  12288 -> w_64x64x96.pmd * B running
###  16384 -> w_64x64x128.pmd *
###  32768 -> w_128x128x64.pmd 

#SBATCH --ntasks=2048

## to be defined !!!!!!!!! 
#SBATCH --time=00:27:43

#SBATCH --job-name="TRIP_PRD_3D"

#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --account=c40
#SBATCH --partition=normal
#SBATCH --ntasks-per-socket=46
#SBATCH --exclusive



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

export MPICH_GPU_SUPPORT_ENABLED=1

ulimit -c 0
ulimit -l unlimited
ulimit -a

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

source ${SCRIPT_DIR}/pfs_evn/lustre_fs_cscs.sh

# ReMap3D diagnostics: enable for one run to validate all-to-all imbalance.
export TRIP_REMAP_MPI_DIAG=1
# Optional timing lines from ReMap3D pack/MPI/unpack sections.
export TRIP_REMAP_TIMING=1

# export MPICH_ALLTOALL_SHORT_MSG_SIZE=0   # tutto trattato come "large"
# export MPICH_ALLTOALL_SYNC_FREQ=16    # prova 8, 16, 32
# export FI_OFI_RXM_SAR_LIMIT=65536    # segment and reassemble limit
# export FI_MR_CACHE_MAX_COUNT=0        # disabilita MR cache (può aiutare a scala)
# export MPICH_OFI_NIC_POLICY=NUMA      # usa NIC più vicina al processo


# Define arguments as a bash array
ARGS=(
    "--yaml_config" "$HOME/git/TRIP/bin/config_scal_acc_alps.yml"
    "--slurm_job_id" "$SLURM_JOB_ID"
    "-ksp_monitor"
    "-ksp_view"
)

echo "Running  with MPS wrapper:"
echo ""
echo ""
echo "Starting TRIP ...... "
echo ""

# Define command as a bash array
SRUN_CMD=(
    "srun"
    # "--switches=1"
    "--cpu-bind=socket"
    "${SCRIPT_DIR}/mps-wrapper.sh"
    "${APP_PATH}/solar_3D"
    "${ARGS[@]}"
)

# Print the command
echo "Command: ${SRUN_CMD[@]}"
echo ""

# Execute the command
"${SRUN_CMD[@]}"
