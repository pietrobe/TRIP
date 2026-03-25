printf "Loading Lustre environment variables...\n"

# =============================================================================
# Filesystem
# =============================================================================
export HDF_FS_TYPE=lustre

# alignment=0 + meta_block_size piccolo — combinazione sicura su Lustre
export HDF_ALIGNMENT_THRESHOLD=0
export HDF_ALIGNMENT_VALUE=4194304     # irrilevante con threshold=0, tenuto per reference
export HDF_META_BLOCK_SIZE=4096        # ← FIX 1: era 4194304, overlap con alignment=0

export HDF_COLL_METADATA_OPS=0
export HDF_COLL_METADATA_WRITE=0
export HDF5_USE_FILE_LOCKING=FALSE

# =============================================================================
# ROMIO collective buffering
# =============================================================================
export HDF_CB_NODES=auto
export HDF_CB_BUFFER_SIZE=134217728
export HDF_ROMIO_CB_WRITE=enable
export HDF_ROMIO_DS_WRITE=disable
export HDF_ROMIO_DS_READ=disable

# Lustre hints
export HDF_LUSTRE_STRIPING_FACTOR=16
export HDF_LUSTRE_STRIPING_UNIT=4194304

# ← FIX 2: GPFS vars removed — not safe to have in environment on Lustre
# export HDF_GPFS_RECV_QUEUES=4          # removed
# export HDF_GPFS_ACCESS_STYLE=write_once  # removed — can cause segfault on Lustre

# =============================================================================
# Cray MPICH — Slingshot (Alps)
# =============================================================================
export MPICH_OFI_STARTUP_CONNECT=1
export MPICH_OFI_NIC_POLICY=BLOCK
export MPICH_OFI_EAGER_THRESHOLD=262144