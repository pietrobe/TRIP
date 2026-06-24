#!/bin/bash
#SBATCH --job-name=hdf5_zslice
#SBATCH --ntasks=12288
#SBATCH --time=00:40:00
#SBATCH --partition=normal
#SBATCH --account=sm111

# =============================================================================
# Output directory — VAST does not use lfs setstripe
# =============================================================================
OUTPUT_DIR=$SCRATCH_NEW/hdf/testc
mkdir -p ${OUTPUT_DIR}
# lfs setstripe -c 32 -S 4m -i -1 ${OUTPUT_DIR}   # Lustre only — not needed on VAST

export HDF_OUTPUT_DIR=${OUTPUT_DIR}
export MPICH_GPU_SUPPORT_ENABLED=0

# =============================================================================
# Force ROMIO to use the UFS driver for VAST
# /ritom/ is the VAST mount prefix on this system
# =============================================================================
export ROMIO_FSTYPE_FORCE='ufs:/ritom/'

# =============================================================================
# ROMIO hints via Cray MPICH — applied to all files opened by this process.
# These are the MPICH-level defaults; the HDF5 FAPL (build_fapl_from_env)
# will override any keys it sets explicitly via MPI_Info.
#
# romio_cb_alltoall : use MPI_Alltoall for aggregation (faster on Slingshot)
# cb_nodes          : one aggregator per node (= SLURM_NNODES)
# romio_ds_write    : disable — VAST writes in-place, no read-before-write
# =============================================================================

export MPIIO_COMMON='romio_cb_read=enable:romio_cb_write=enable:romio_cb_alltoall=enable'
# romio_ds_write/read are left off here; the HDF5 FAPL (build_fapl_from_env)
# explicitly sets romio_ds_write=disable / romio_ds_read=disable for VAST.
# Enabling them here would conflict and the interplay is MPICH-implementation-specific.
# cb_config_list removed: it overrides cb_nodes and forces all tasks as aggregators
# regardless of HDF_CB_NODES. Let the FAPL cb_nodes setting take effect instead.

export MPICH_MPIIO_HINTS="*:${MPIIO_COMMON}"

# =============================================================================
# Problem dimensions
# =============================================================================
export HDF_N_X=128
export HDF_N_Y=128
export HDF_N_Z=292

export HDF_N_FREQUENCIES=96
export HDF_N_INCL=8
export HDF_N_AZIMUTH=16
export HDF_N_STOKES=4
export HDF_NORMALIZE_OUTPUT=0

# step_z=16: N_z/step_z = 143/16 ≈ 8.9375 H5Dwrite_multi collective rounds (was 32 with step_z=4).
# Each round costs a Slingshot alltoall + VAST write; 4× fewer rounds = 4× less overhead.
# chunk_z MUST equal step_z — writes land on exact chunk boundaries (no partial-chunk I/O).
# chunk = {1, 16, 16, 16, 8, 96} × float32 = 12 MB — large sequential blocks for VAST NVMe.
# Memory per rank: 4 Stokes × (Nlx × Nly × 16 × 16 × 8 × 96 × 4B) ≈ 2.3 GB — fine on Alps.
export HDF_Z_SLICE_WIDTH=4
export HDF_STEP_Z=4
export HDF_CHUNK_Z=4   # must equal HDF_STEP_Z

# =============================================================================
# HDF5 / ROMIO tuning — VAST Universal Storage
# =============================================================================
export HDF_FS_TYPE=vast

# ROMIO collective buffering — all tasks as aggregators for VAST.
# With cb_nodes=auto (≈32 nodes), 256 ranks must first redistribute data
# to 32 aggregators via alltoall, then aggregators write to VAST.
# With cb_nodes=NTASKS=256, every rank is its own aggregator: the alltoall
# becomes a no-op (send-to-self), and all 256 ranks write to VAST
# simultaneously — leveraging VAST's native high-concurrency NVMe backend.
export HDF_CB_NODES=${SLURM_NNODES}    # 256 → direct writes, no alltoall redistribution
export HDF_CB_BUFFER_SIZE=536870912    # 512 MB — enough for 192 MB/Stokes per rank at step_z=16
export HDF_ROMIO_CB_WRITE=enable
export HDF_ROMIO_DS_WRITE=disable      # no read-before-write on VAST
export HDF_ROMIO_DS_READ=disable

# HDF5 object alignment — fully disabled on VAST.
#
# Natural I/O block : 96 × step_z × N_y × N_az × sizeof(float32)
#                   = 96 × 16 × 128 × 8 × 4 = 6,291,456 bytes = 6 MB
# HDF5 chunk        : {1, 16, step_z, 16, 8, 96} × sizeof(float32)
#                   = 1 × 16 × 16 × 16 × 8 × 96 × 4 = 12,582,912 bytes = 12 MB (= 2 natural blocks)
#
# Any alignment value that does NOT divide 3 MB (e.g., 4 MB) pads every chunk
# and inflates the file proportionally.  Setting threshold above the largest
# possible object (1 TB) makes alignment a guaranteed no-op regardless of
# the value setting.
export HDF_ALIGNMENT_THRESHOLD=1099511627776  # 1 TB — no object ever qualifies; alignment is disabled
export HDF_ALIGNMENT_VALUE=1                  # safety net: if anything qualified, pad to 1 byte = 0 waste

# Metadata block size = 1 natural block (1.5 MB).
# Keeps metadata aggregated (fewer small writes) while staying commensurate
# with the data block granularity and avoiding waste at the end of the last block.
export HDF_META_BLOCK_SIZE=1572864            # 1.5 MB = natural block
export HDF_COLL_METADATA_OPS=1
export HDF_COLL_METADATA_WRITE=1
export HDF_META_READ_ATTEMPTS=20       # default is 1; retries on transient checksum failures

export HDF5_USE_FILE_LOCKING=FALSE

# =============================================================================
# Lustre settings — commented out, not applicable on VAST
# =============================================================================
# export HDF_FS_TYPE=lustre
# export HDF_ALIGNMENT_THRESHOLD=0       # 0 = disable alignment on Lustre
# export HDF_ALIGNMENT_VALUE=4194304     # 4 MB stripe unit
# export HDF_META_BLOCK_SIZE=4194304
# export HDF_COLL_METADATA_OPS=0         # avoid deadlock with slab_comm on Cray MPICH
# export HDF_COLL_METADATA_WRITE=0
# export HDF_CB_BUFFER_SIZE=134217728    # 128 MB
# export HDF_LUSTRE_STRIPING_FACTOR=16   # must match lfs setstripe -c
# export HDF_LUSTRE_STRIPING_UNIT=4194304
# export HDF_GPFS_RECV_QUEUES=4
# export HDF_GPFS_ACCESS_STYLE=write_once

# =============================================================================
# Cray MPICH — Slingshot (Alps)
# =============================================================================
export MPICH_OFI_STARTUP_CONNECT=1
export MPICH_OFI_NIC_POLICY=BLOCK
export MPICH_OFI_EAGER_THRESHOLD=262144   # 256 KB eager limit

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
    printf \"File size:  %.2f GB  (4 datasets × %.2f GB each)\n\", 4*one, one
}"
echo "================================================================"

# =============================================================================
# Run
# =============================================================================
srun -A c40 ./hdf_output_mpi_example

# =============================================================================
# Post-run
# =============================================================================
echo "Output files:"
ls -lh ${OUTPUT_DIR}/*.h5 2>/dev/null || echo "  (no .h5 files found)"
