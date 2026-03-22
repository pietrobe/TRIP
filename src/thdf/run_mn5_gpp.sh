#!/bin/bash
#SBATCH --job-name=hdf5_zslice
#SBATCH --partition=gpp
#SBATCH --qos=gp_ehpc
#SBATCH --account=ehpc238
### # SBATCH --nodes=447
### # SBATCH --ntasks-per-node=112
### # SBATCH --cpus-per-task=1
#SBATCH --ntasks=1024
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%x_%j.out
#SBATCH --error=slurm-%x_%j.err
#SBATCH --exclusive

export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
ulimit -s unlimited

BSC_GROUP=ehpc238
OUTPUT_DIR=/gpfs/scratch/${BSC_GROUP}/hdf5_output/
mkdir -p ${OUTPUT_DIR}
export HDF_OUTPUT_DIR=${OUTPUT_DIR}

export HDF_N_X=64
export HDF_N_Y=64
export HDF_N_Z=64
export HDF_N_FREQUENCIES=96
export HDF_N_INCL=16
export HDF_N_AZIMUTH=8
export HDF_N_STOKES=4
export HDF_NORMALIZE_OUTPUT=0
export HDF_Z_SLICE_WIDTH=2

export HDF_FS_TYPE=gpfs
# export HDF_META_BLOCK_SIZE=8388608
export HDF_COLL_METADATA_OPS=0
export HDF_COLL_METADATA_WRITE=0
export HDF5_USE_FILE_LOCKING=FALSE

# ✅ Soglia 1 MB — allinea solo i chunk (> 384 KB), non i metadata piccoli
# export HDF_ALIGNMENT_THRESHOLD=1048576   # 1 MB
# export HDF_ALIGNMENT_VALUE=4194304       # 4 MB (laptop) | 8388608 (MN5 GPFS)

# ✅ meta_block_size piccolo — evita il gap di 165 MB
# export HDF_META_BLOCK_SIZE=4096          # 4 KB — nessun pre-allocazione

export HDF_ALIGNMENT_THRESHOLD=0     # disabilita alignment — elimina unaccounted space
export HDF_ALIGNMENT_VALUE=1         # irrilevante con threshold=0
export HDF_META_BLOCK_SIZE=8388608   # 8 MB — necessario su GPFS, causa lo stallo se rimosso

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

n_slices=$(( (HDF_N_Z + HDF_Z_SLICE_WIDTH - 1) / HDF_Z_SLICE_WIDTH ))
echo "================================================================"
echo "Job:        ${SLURM_JOB_ID}"
echo "Nodes:      ${SLURM_JOB_NUM_NODES}"
echo "Tasks:      ${SLURM_NTASKS}"
echo "Output:     ${OUTPUT_DIR}"
echo "Z slices:   ${n_slices} x ${HDF_Z_SLICE_WIDTH} planes"
awk "BEGIN {printf \"File size:  %.2f GB/slice\n\", 128 * 128 * 16 * 16 * 8 * 96 * 4 * 4 / 1e9}"
echo "Total size: ~384 GB (4 datasets x 96 GB)"
echo "================================================================"

srun --ntasks=${SLURM_NTASKS}    \
     --ntasks-per-node=112       \
     --cpus-per-task=1           \
     ./hdf_output_mpi_example

echo "Output files:"
ls -lh ${OUTPUT_DIR}/3d_field/*.h5 2>/dev/null || echo "  (no .h5 files found)"
mmlsquota -u ${USER} --block-size auto gpfs_scratch 2>/dev/null || true