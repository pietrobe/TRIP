#!/bin/bash
# =============================================================================
# GPFS_fs.sh — IBM Spectrum Scale (GPFS) tuning for HDF5 parallel I/O
#
# Source this file before launching the application:
#   source GPFS_fs.sh
#
# Tested on:
#   MareNostrum 5 (BSC) — Intel MPI 2021, GPFS /gpfs/scratch
#
# All variables can be overridden after sourcing by re-exporting them.
# =============================================================================

# =============================================================================
# ROMIO hints (Cray MPICH: MPICH_MPIIO_HINTS / Intel MPI: ROMIO_HINTS file)
#
# romio_cb_read/write : enable collective buffering
# romio_ds_write      : DISABLED — avoids read-modify-write overhead on GPFS
# romio_ds_read       : DISABLED
#
# NOTE: cb_config_list is intentionally NOT set here.
#   Setting cb_config_list=#*:*# overrides cb_nodes and forces all tasks as
#   aggregators, eliminating collective buffering benefits on GPFS.
# =============================================================================
_MPIIO_HINTS='romio_cb_read=enable'
_MPIIO_HINTS="${_MPIIO_HINTS}:romio_cb_write=enable"
_MPIIO_HINTS="${_MPIIO_HINTS}:romio_ds_write=disable"
_MPIIO_HINTS="${_MPIIO_HINTS}:romio_ds_read=disable"

export MPICH_MPIIO_HINTS="*:${_MPIIO_HINTS}"
unset _MPIIO_HINTS

# =============================================================================
# HDF5 filesystem type
# Selects the GPFS branch in build_fapl_from_env().
# =============================================================================
export HDF_FS_TYPE=gpfs

# =============================================================================
# ROMIO collective buffering — one aggregator per physical node
#
# cb_nodes=auto resolves to the physical node count via MPI_Comm_split_type.
# One aggregator per node aligns with the GPFS NSD server topology and avoids
# over-subscribing the GPFS lock manager.
# =============================================================================
export HDF_CB_NODES=auto
export HDF_CB_BUFFER_SIZE=268435456    # 256 MB per aggregator

export HDF_ROMIO_CB_WRITE=enable
export HDF_ROMIO_DS_WRITE=disable
export HDF_ROMIO_DS_READ=disable

# =============================================================================
# GPFS-specific ROMIO hints (passed via MPI_Info in build_fapl_from_env)
#
# HDF_GPFS_RECV_QUEUES  : NOT SET — causes MPI_Abort on MN5 GPFS client.
# HDF_GPFS_ACCESS_STYLE : NOT SET — "write_once" conflicts with the serial
#   H5Fopen that rank 0 performs for metadata (Phase 2); ROMIO caches regions
#   as write-only and rank 0's re-read gets stale data.
# HDF_GPFS_NUM_IO_NODES : optional — 0 (default) lets GPFS choose.
# =============================================================================

# =============================================================================
# HDF5 object alignment
#
# Disabled (threshold=0) to avoid unaccounted space in the file.
# GPFS block sizes range from 512 KB to 4 MB; aligning HDF5 objects to that
# boundary pads every chunk.  Enable only if profiling shows a benefit:
#   export HDF_ALIGNMENT_THRESHOLD=1048576   # 1 MB
#   export HDF_ALIGNMENT_VALUE=4194304       # 4 MB (match GPFS block size)
# =============================================================================
export HDF_ALIGNMENT_THRESHOLD=0      # 0 = disable alignment entirely
export HDF_ALIGNMENT_VALUE=1          # irrelevant when threshold=0

# =============================================================================
# HDF5 metadata block size
#
# 8 MB is required on GPFS to prevent the HDF5 metadata aggregator from
# issuing many small writes that stall the GPFS client.  Do not reduce
# below 4 MB; the library may deadlock on metadata flush on some GPFS builds.
# =============================================================================
export HDF_META_BLOCK_SIZE=8388608    # 8 MB

# =============================================================================
# HDF5 collective metadata I/O — DISABLED on GPFS
#
# With HDF_COLL_METADATA_OPS=1, HDF5 issues collective MPI calls for every
# metadata read/write.  This deadlocks on Cray MPICH and Intel MPI when a
# sub-communicator (e.g. the slab_comm used in Phase 2) opens the file while
# the remaining ranks are blocked in MPI_Barrier on MPI_COMM_WORLD.
# GPFS provides its own coherency via the distributed lock manager (DLM), so
# collective HDF5 metadata is not needed for correctness.
# =============================================================================
export HDF_COLL_METADATA_OPS=0
export HDF_COLL_METADATA_WRITE=0

# GPFS DLM ensures coherency; a low retry count is sufficient.
export HDF_META_READ_ATTEMPTS=5

# =============================================================================
# HDF5 file locking
#
# Disabled: MPI collectives coordinate parallel access within the job.
# Disabling avoids contention on the GPFS DLM when thousands of ranks open
# the same file simultaneously.
# =============================================================================
export HDF5_USE_FILE_LOCKING=FALSE
