printf "Loading GPFS environment variables for MN5...\n"


export HDF_FS_TYPE=gpfs
# export HDF_META_BLOCK_SIZE=8388608
export HDF_COLL_METADATA_OPS=0
export HDF_COLL_METADATA_WRITE=0
export HDF5_USE_FILE_LOCKING=FALSE

# Threshold 1 MB — align only the chunk (> 384 KB), not the small metadata
# export HDF_ALIGNMENT_THRESHOLD=1048576   # 1 MB
# export HDF_ALIGNMENT_VALUE=4194304       # 4 MB (laptop) | 8388608 (MN5 GPFS)

# meta_block_size small — prevent the gap of 165 MB
# export HDF_META_BLOCK_SIZE=4096          # 4 KB — no pre-allocation

export HDF_ALIGNMENT_THRESHOLD=0     # disable alignment — eliminate unaccounted space
export HDF_ALIGNMENT_VALUE=1         # irrelevant with threshold=0
export HDF_META_BLOCK_SIZE=8388608   # 8 MB — necessary on GPFS, causes stall if removed

export HDF_CB_NODES=auto
export HDF_CB_BUFFER_SIZE=268435456
export HDF_ROMIO_CB_WRITE=enable
export HDF_ROMIO_DS_WRITE=disable
export HDF_ROMIO_DS_READ=disable
# HDF_GPFS_RECV_QUEUES: removed — causes MPI abort on MN5 GPFS client
# HDF_GPFS_ACCESS_STYLE: removed — conflicts with serial H5Fopen for metadata

export I_MPI_OFI_PROVIDER=mlx
export I_MPI_FABRICS=shm:ofi
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core
export I_MPI_ADJUST_ALLREDUCE=4
export I_MPI_ADJUST_BARRIER=4
export I_MPI_ADJUST_BCAST=4
export I_MPI_OFI_STARTUP_CONNECT=1

ROMIO_HINTS_FILE=${OUTPUT_DIR}/romio_hints.txt
cat > ${ROMIO_HINTS_FILE} << EOF
romio_cb_write      enable
romio_ds_write      disable
romio_ds_read       disable
cb_buffer_size      268435456
EOF
export ROMIO_HINTS=${ROMIO_HINTS_FILE}

export EAR_POWER_POLICY=monitoring