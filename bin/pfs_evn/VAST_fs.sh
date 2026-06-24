#!/bin/bash
# =============================================================================
# VAST_fs.sh — VAST Universal Storage tuning for Cray MPICH / HDF5 on Alps
#
# Source this file before launching the application:
#   source VAST_fs.sh
#
# All variables can be overridden after sourcing by re-exporting them.
# =============================================================================

# =============================================================================
# Force ROMIO to use the UFS driver for VAST
# VAST exposes a POSIX-compatible interface; the UFS driver is the correct
# ROMIO backend for the /ritom/ mount point on this system.
# =============================================================================
export ROMIO_FSTYPE_FORCE='ufs:/ritom/'

# =============================================================================
# ROMIO hints (Cray MPICH MPICH_MPIIO_HINTS)
#
# romio_cb_read/write : enable collective buffering for reads and writes
# romio_cb_alltoall   : use MPI_Alltoall for aggregation (faster on Slingshot)
# romio_ds_write      : DISABLED — VAST writes in-place; no read-before-write
# romio_ds_read       : DISABLED — avoid data sieving on VAST
#
# NOTE: cb_config_list is intentionally NOT set here.
#   Setting cb_config_list=#*:*# forces all ranks as aggregators and silently
#   overrides any cb_nodes value set by the HDF5 FAPL or HDF_CB_NODES below.
#   Omitting it allows the FAPL cb_nodes (= SLURM_NNODES) to take effect.
# =============================================================================
_MPIIO_HINTS='romio_cb_read=enable'
_MPIIO_HINTS="${_MPIIO_HINTS}:romio_cb_write=enable"
_MPIIO_HINTS="${_MPIIO_HINTS}:romio_cb_alltoall=enable"
_MPIIO_HINTS="${_MPIIO_HINTS}:romio_ds_write=disable"
_MPIIO_HINTS="${_MPIIO_HINTS}:romio_ds_read=disable"

export MPICH_MPIIO_HINTS="*:${_MPIIO_HINTS}"
unset _MPIIO_HINTS

# =============================================================================
# HDF5 filesystem type
# Selects the Lustre/GPFS/VAST branch in build_fapl_from_env().
# =============================================================================
export HDF_FS_TYPE=vast

# =============================================================================
# ROMIO collective buffering — cb_nodes = one aggregator per physical node
#
# cb_nodes=auto resolves to the number of physical nodes in the communicator
# (detected via MPI_Comm_split_type in build_fapl_from_env).  With one
# aggregator per node, each aggregator writes a large, node-sized stripe to
# VAST, maximising per-aggregator bandwidth without over-subscribing VAST's
# metadata server.
#
# Alternative: set to a fixed number, e.g. export HDF_CB_NODES=64
# =============================================================================
export HDF_CB_NODES=auto

# cb_buffer_size per aggregator: 512 MB is adequate for VAST's NVMe backend.
# Increase if cb_nodes is small and the per-aggregator write stripe is large.
export HDF_CB_BUFFER_SIZE=536870912    # 512 MB

export HDF_ROMIO_CB_WRITE=enable
export HDF_ROMIO_DS_WRITE=disable
export HDF_ROMIO_DS_READ=disable

# =============================================================================
# HDF5 object alignment
#
# VAST places data internally on its NVMe object store; there is no fixed
# stripe unit to align to.  Setting threshold above any realistic object size
# disables alignment padding, avoiding file inflation.
# =============================================================================
export HDF_ALIGNMENT_THRESHOLD=1099511627776  # 1 TB — effectively disables alignment
export HDF_ALIGNMENT_VALUE=1                  # 1 byte = no padding if anything qualifies

# =============================================================================
# HDF5 metadata block size
#
# Aggregates small metadata writes into 1.5 MB blocks.  Commensurate with a
# typical z-slab I/O block size; reduces metadata round-trips to VAST.
# Tune alongside HDF_STEP_Z if the natural I/O block size changes.
# =============================================================================
export HDF_META_BLOCK_SIZE=1572864            # 1.5 MB

# =============================================================================
# HDF5 collective metadata I/O
#
# Both ops and write are enabled so that all metadata reads and writes go
# through a single MPI collective (rank-0 broadcasts), preventing races on
# the fixed-array chunk index pages under high parallelism.
# =============================================================================
export HDF_COLL_METADATA_OPS=1
export HDF_COLL_METADATA_WRITE=1

# Retry count for metadata page reads before declaring a checksum failure.
# Default is 1 (non-SWMR); 20 absorbs transient VAST client-cache glitches.
export HDF_META_READ_ATTEMPTS=20

# =============================================================================
# HDF5 file locking
#
# Disabled for parallel HDF5: coordination is via MPI collectives, not POSIX
# locks.  POSIX locks on VAST can cause spurious failures with many openers.
# =============================================================================
export HDF5_USE_FILE_LOCKING=FALSE

# =============================================================================
# Cray MPICH / Slingshot (Alps) settings
#
# MPICH_OFI_STARTUP_CONNECT : pre-connect all OFI endpoints at MPI_Init;
#   avoids connection latency during the first alltoall of the write phase.
# MPICH_OFI_NIC_POLICY=BLOCK : assign NICs in contiguous blocks per node,
#   reducing cross-node NIC sharing and improving injection bandwidth.
# MPICH_OFI_EAGER_THRESHOLD  : messages below this size use the eager
#   protocol (no rendezvous); 256 KB matches the typical aggregator stripe.
# =============================================================================
export MPICH_OFI_STARTUP_CONNECT=1
export MPICH_OFI_NIC_POLICY=BLOCK
export MPICH_OFI_EAGER_THRESHOLD=262144       # 256 KB
