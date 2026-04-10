#!/bin/bash
#SBATCH --job-name=hdf5_zslice
#SBATCH --ntasks=192 
#SBATCH --time=00:20:00
#SBATCH --partition=debug
#SBATCH --account=c40

# =============================================================================
# Lustre stripe — must be set BEFORE the files are created
# =============================================================================
OUTPUT_DIR=$SCRATCH/hdf_test
mkdir -p ${OUTPUT_DIR}
lfs setstripe -c 16 -S 4m -i -1 ${OUTPUT_DIR}

export HDF_OUTPUT_DIR=${OUTPUT_DIR}
export MPICH_GPU_SUPPORT_ENABLED=0

# =============================================================================
# Problem dimensions
# =============================================================================
export HDF_N_X=32
export HDF_N_Y=32
export HDF_N_Z=32
export HDF_N_FREQUENCIES=96
export HDF_N_INCL=16
export HDF_N_AZIMUTH=8
export HDF_N_STOKES=4
export HDF_NORMALIZE_OUTPUT=0

# HDF_STEP_Z controls the write buffer size per iteration — different parameter.
# Both must be set. For single-file mode also set HDF_CHUNK_Z=HDF_STEP_Z.
export HDF_Z_SLICE_WIDTH=4
export HDF_STEP_Z=4            # write buffer: one full slice per iteration
export HDF_CHUNK_Z=4           # must equal HDF_STEP_Z for single-file mode

# =============================================================================
# Filesystem
# =============================================================================
export HDF_FS_TYPE=lustre

# HDF_ALIGNMENT_THRESHOLD=1 causes ~10 GB unaccounted space per file
# (every object including tiny metadata headers is padded to 4 MB).
# Set threshold to 0 to disable alignment on Lustre with chunked HDF5 —
# ROMIO collective buffering already handles alignment internally.
export HDF_ALIGNMENT_THRESHOLD=0      # was 1 → caused 13x file size inflation
export HDF_ALIGNMENT_VALUE=4194304    # 4 MB stripe unit (kept for reference)
export HDF_META_BLOCK_SIZE=4194304    # 4 MB metadata buffer — keep on Lustre

# HDF_COLL_METADATA_OPS=1 causes deadlock on subcommunicators
# with Cray MPICH on Alps when slab_comm != MPI_COMM_WORLD.
export HDF_COLL_METADATA_OPS=0        # was 1 → deadlock with slab_comm
export HDF_COLL_METADATA_WRITE=0      # was 1 → idem

export HDF5_USE_FILE_LOCKING=FALSE

# =============================================================================
# ROMIO collective buffering
# =============================================================================
export HDF_CB_NODES=auto
export HDF_CB_BUFFER_SIZE=134217728    # 128 MB per aggregator
export HDF_ROMIO_CB_WRITE=enable
export HDF_ROMIO_DS_WRITE=disable
export HDF_ROMIO_DS_READ=disable

# Lustre-specific hints
export HDF_LUSTRE_STRIPING_FACTOR=16
export HDF_LUSTRE_STRIPING_UNIT=4194304    # 4 MB — must match lfs setstripe -S

# GPFS hints (unused on Lustre — kept for portability)
export HDF_GPFS_RECV_QUEUES=4
export HDF_GPFS_ACCESS_STYLE=write_once

# =============================================================================
# Cray MPICH — Slingshot (Alps)
# =============================================================================
export MPICH_OFI_STARTUP_CONNECT=1
export MPICH_OFI_NIC_POLICY=BLOCK     # NUMA removed: requires per-NUMA confinement

# UCX_TLS is irrelevant on Alps — Cray MPICH uses OFI/CXI directly,
# not UCX. These variables are silently ignored but can cause confusion.
# Remove UCX_TLS, UCX_RC_TIMEOUT, UCX_UD_TIMEOUT.
# Use MPICH_OFI_* variables instead for timeout tuning.
export MPICH_OFI_EAGER_THRESHOLD=262144   # 256 KB eager limit for large jobs

# Generous timeout for H5Fclose on large files (384 GB flush can take >60s)
# Cray MPICH does not have a direct timeout env var — increase job time instead.
# The --time=00:20:00 above is likely too short for 24576 ranks + 384 GB.
# Recommendation: increase to at least 01:00:00 for production runs.

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
awk "BEGIN {printf \"File size:  %.2f GB/slice\n\", ${HDF_N_X} * ${HDF_N_Y} * ${HDF_Z_SLICE_WIDTH} * ${HDF_N_INCL} * ${HDF_N_AZIMUTH} * ${HDF_N_FREQUENCIES} * 4 * 4 / 1e9}"
echo "================================================================"

# =============================================================================
# Run
# =============================================================================
srun ./hdf_output_mpi_example

# =============================================================================
# Post-run
# =============================================================================
echo "Output files:"
ls -lh ${OUTPUT_DIR}/*.h5 2>/dev/null || echo "  (no .h5 files found)"
