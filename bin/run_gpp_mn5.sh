#!/bin/bash -l
#SBATCH --job-name="TRIP_PRD_3D"
#SBATCH --account=ehpc597
#SBATCH --qos=gp_ehpc
#SBATCH --time=05:47:28
#### #SBATCH --nodes=476
#SBATCH --ntasks=32768
#### #SBATCH --ntasks=
#SBATCH --ntasks-per-node=112
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --propagate=STACK,MEMLOCK
#### SBATCH --exclude=gs21r1b06


export APP_PATH=${HOME}/git/TRIP/build
export SCRIPT_DIR=${HOME}/git/TRIP/bin

echo "=========================================="
echo "User-Defined Variables:"
echo "=========================================="
printf "%-30s %s\n" "APP_PATH:"   "$APP_PATH"
printf "%-30s %s\n" "SCRIPT_DIR:" "$SCRIPT_DIR"
echo "=========================================="
echo ""

echo "Starting TRIP ...... "
echo ""

source ${SCRIPT_DIR}/pfs_evn/GPFS_mn5.sh
source ${SCRIPT_DIR}/slurm_env_print.sh

export RII_AA_SCRATCH_DIR=$TMPDIR
export RII_AA_JOB_PATTERN="TRIP_job_$SLURM_JOB_ID"
echo "RII_AA_SCRATCH_DIR: $RII_AA_SCRATCH_DIR"
echo "RII_AA_JOB_PATTERN: $RII_AA_JOB_PATTERN"
echo ""

# ----------------------------------------------------------------------
# UCX / MPI tuning for 65k-rank scale
# ----------------------------------------------------------------------
# Use Dynamically Connected transport — RC-reliable, scales past UD's
# wireup ceiling, avoids UD timeout cascades at 64k+ ranks.
export UCX_TLS=dc_x,sm,self

# Memory registration: rcache (default, but explicit) + skip CUDA
# memtype hooks (CPU-only job) + cheaper rcache lookups.
export UCX_IB_REG_METHODS=rcache
export UCX_MEMTYPE_CACHE=n
export UCX_RCACHE_CHECK_PFN=0
# At 112 ranks/node the per-process rcache accumulates enough pinned VMAs to
# exhaust vm.max_map_count (~65530), causing ibv_reg_mr ENOMEM on tiny buffers.
# Huge pages collapse ~512 4K-page VMAs into one 2M-page VMA each (primary fix).
# The region cap provides a hard backstop if huge pages aren't available.
export UCX_IB_HUGETLB=y
export UCX_RCACHE_MAX_REGIONS=2048

# Intel MPI: do NOT enable async progress threads with 112 ranks/node
# and cpus-per-task=1 — they contend with the ranks themselves.
unset I_MPI_ASYNC_PROGRESS

# Rabenseifner allreduce, good for medium/large messages at scale.
export I_MPI_ADJUST_ALLREDUCE=5



ARGS=(
    "--yaml_config" "$HOME/git/TRIP/bin/config_mn5.yml"
    "-ksp_monitor"
    "-ksp_view"
)

echo "Executing command:"
echo " srun  ${APP_PATH}/solar_3D ${ARGS[*]}"
echo ""

srun ${APP_PATH}/solar_3D "${ARGS[@]}"