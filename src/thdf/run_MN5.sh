#!/bin/bash
#SBATCH --job-name=hdf5_zslice
#SBATCH --partition=gpp
#SBATCH --qos=<QOS>
#SBATCH --account=<ACCOUNT>
#SBATCH --ntasks=<NTASKS>
#SBATCH --ntasks-per-node=112
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%x_%j.out
#SBATCH --error=slurm-%x_%j.err
#SBATCH --exclusive

ulimit -s unlimited
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

# =============================================================================
# Output directory — GPFS scratch
# =============================================================================
OUTPUT_DIR=/gpfs/scratch/<GROUP>/<USER>/hdf/test
mkdir -p ${OUTPUT_DIR}
export HDF_OUTPUT_DIR=${OUTPUT_DIR}

# =============================================================================
# GPFS / ROMIO / HDF5 environment
# =============================================================================
source "$(dirname "$0")/GPFS_fs.sh"

# =============================================================================
# Intel MPI / InfiniBand (MN5 — Sapphire Rapids, HDR InfiniBand)
# =============================================================================
export I_MPI_OFI_PROVIDER=mlx
export I_MPI_FABRICS=shm:ofi
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core
export I_MPI_OFI_STARTUP_CONNECT=1

# =============================================================================
# EAR energy monitoring (MN5)
# =============================================================================
export EAR_POWER_POLICY=monitoring

# =============================================================================
# Problem dimensions
# =============================================================================
export HDF_N_X=128
export HDF_N_Y=128
export HDF_N_Z=256

export HDF_N_FREQUENCIES=96
export HDF_N_INCL=16
export HDF_N_AZIMUTH=8
export HDF_N_STOKES=4
export HDF_NORMALIZE_OUTPUT=0

export HDF_Z_SLICE_WIDTH=4
export HDF_STEP_Z=4

# =============================================================================
# Print job info
# =============================================================================
n_slices=$(( (HDF_N_Z + HDF_Z_SLICE_WIDTH - 1) / HDF_Z_SLICE_WIDTH ))
echo "================================================================"
echo "Job:        ${SLURM_JOB_ID}"
echo "Nodes:      ${SLURM_JOB_NUM_NODES}"
echo "Tasks:      ${SLURM_NTASKS}"
echo "Output:     ${OUTPUT_DIR}"
echo "Z slices:   ${n_slices} x ${HDF_Z_SLICE_WIDTH} planes"
awk "BEGIN {
    one = ${HDF_N_X} * ${HDF_N_Y} * ${HDF_N_Z} * ${HDF_N_INCL} * ${HDF_N_AZIMUTH} * ${HDF_N_FREQUENCIES} * 4 / 1e9
    printf \"File size:  %.2f GB  (4 datasets x %.2f GB each)\n\", 4*one, one
}"
echo "================================================================"

# =============================================================================
# Run
# =============================================================================
srun --ntasks=${SLURM_NTASKS} \
     --ntasks-per-node=112    \
     --cpus-per-task=1        \
     ./hdf_output_mpi_example

# =============================================================================
# Post-run
# =============================================================================
echo "Output files:"
ls -lh ${OUTPUT_DIR}/*.h5 2>/dev/null || echo "  (no .h5 files found)"
mmlsquota --block-size auto gpfs_scratch 2>/dev/null || true
