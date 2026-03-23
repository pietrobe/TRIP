
printf "Loading Lustre environment variables...\n"

# OUTPUT_DIR=$SCRATCH/hdf_test
# mkdir -p ${OUTPUT_DIR}
# lfs setstripe -c 16 -S 4m -i -1 ${OUTPUT_DIR}

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