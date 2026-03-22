/**
 * hdf5_fapl_from_env.c
 *
 * Builds an HDF5 File Access Property List (FAPL) entirely from
 * environment variables. No recompilation needed to tune for
 * Lustre, GPFS, or different job sizes.
 *
 * Usage in job script:
 *
 *   # Required
 *   export HDF_FS_TYPE=lustre          # or gpfs
 *
 *   # ROMIO collective buffering
 *   export HDF_CB_NODES=auto           # "auto" = 1 per physical node
 *   export HDF_CB_BUFFER_SIZE=134217728
 *   export HDF_ROMIO_CB_WRITE=enable
 *   export HDF_ROMIO_DS_WRITE=disable
 *
 *   # Lustre specific
 *   export HDF_LUSTRE_STRIPING_FACTOR=16
 *   export HDF_LUSTRE_STRIPING_UNIT=4194304
 *
 *   # GPFS specific
 *   export HDF_GPFS_RECV_QUEUES=4
 *   export HDF_GPFS_ACCESS_STYLE=write_once
 *
 *   # HDF5 layer
 *   export HDF_ALIGNMENT_THRESHOLD=1
 *   export HDF_ALIGNMENT_VALUE=4194304
 *   export HDF_META_BLOCK_SIZE=4194304
 *   export HDF_COLL_METADATA_OPS=1
 *   export HDF_COLL_METADATA_WRITE=1
 */

#include "hdf_fapl_from_env.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* =========================================================================
 * Environment helpers
 * ========================================================================= */

/**
 * Read a string env variable.
 * Returns def if not set or empty.
 */
static const char *
env_str(const char *name, const char *def) {
  const char *v = getenv(name);
  return (v && v[0]) ? v : def;
}

/**
 * Read a long integer env variable.
 * Returns def if not set, empty, or unparseable.
 */
static long
env_long(const char *name, long def) {
  const char *v = getenv(name);
  if (!v || !v[0])
    return def;
  char *end;
  long  val = strtol(v, &end, 10);
  return (end != v) ? val : def;
}

/**
 * Read a size_t env variable (used for byte counts).
 */
static size_t
env_size(const char *name, size_t def) {
  const char *v = getenv(name);
  if (!v || !v[0])
    return def;
  char              *end;
  unsigned long long val = strtoull(v, &end, 10);
  return (end != v) ? (size_t)val : def;
}

/**
 * Count physical nodes in a communicator.
 * Uses MPI shared-memory split — portable across any topology.
 */
static int
count_nodes(MPI_Comm comm) {
  MPI_Comm node_comm;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
  int node_rank;
  MPI_Comm_rank(node_comm, &node_rank);
  MPI_Comm_free(&node_comm);

  int is_root = (node_rank == 0) ? 1 : 0;
  int n_nodes;
  MPI_Allreduce(&is_root, &n_nodes, 1, MPI_INT, MPI_SUM, comm);
  return n_nodes;
}

/* =========================================================================
 * build_fapl_from_env
 *
 * All tuning parameters come from environment variables.
 * Falls back to sensible defaults when variables are not set.
 * ========================================================================= */
hid_t
build_fapl_from_env(MPI_Comm comm) {

  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  /* ----------------------------------------------------------------
   * Filesystem type
   * ---------------------------------------------------------------- */
  const char *fs_type = env_str("HDF_FS_TYPE", "lustre");

  /* ----------------------------------------------------------------
   * cb_nodes
   *
   * HDF_CB_NODES=auto  → count physical nodes in comm (recommended)
   * HDF_CB_NODES=N     → use exactly N aggregators
   * ---------------------------------------------------------------- */
  const char *cb_nodes_env = env_str("HDF_CB_NODES", "auto");
  char        cb_nodes_str[32];

  if (strcmp(cb_nodes_env, "auto") == 0) {
    int n_nodes = count_nodes(comm);
    snprintf(cb_nodes_str, sizeof(cb_nodes_str), "%d", n_nodes);
    if (mpi_rank == 0)
      printf("[fapl] HDF_CB_NODES=auto → cb_nodes=%d (physical nodes)\n", n_nodes);
  } else {
    snprintf(cb_nodes_str, sizeof(cb_nodes_str), "%s", cb_nodes_env);
    if (mpi_rank == 0)
      printf("[fapl] HDF_CB_NODES=%s\n", cb_nodes_str);
  }

  /* ----------------------------------------------------------------
   * ROMIO hints
   * ---------------------------------------------------------------- */
  size_t      cb_buffer_size = env_size("HDF_CB_BUFFER_SIZE", 134217728); /* 128 MB */
  const char *romio_cb_write = env_str("HDF_ROMIO_CB_WRITE", "enable");
  const char *romio_ds_write = env_str("HDF_ROMIO_DS_WRITE", "disable");
  const char *romio_ds_read  = env_str("HDF_ROMIO_DS_READ", "disable");

  /* ----------------------------------------------------------------
   * Lustre hints
   * ---------------------------------------------------------------- */
  long   lustre_striping_factor = env_long("HDF_LUSTRE_STRIPING_FACTOR", 16);
  size_t lustre_striping_unit   = env_size("HDF_LUSTRE_STRIPING_UNIT", 4194304); /* 4 MB */

  /* ----------------------------------------------------------------
   * GPFS hints
   * ---------------------------------------------------------------- */
  long        gpfs_recv_queues  = env_long("HDF_GPFS_RECV_QUEUES", 4);
  const char *gpfs_access_style = env_str("HDF_GPFS_ACCESS_STYLE", "write_once");
  long        gpfs_num_io_nodes = env_long("HDF_GPFS_NUM_IO_NODES", 0); /* 0 = not set */

  /* ----------------------------------------------------------------
   * HDF5 layer
   * ---------------------------------------------------------------- */
  // ✅ Allinea solo oggetti >= 1 MB — i chunk vengono allineati,
  //    i metadata piccoli no
  hsize_t alignment_threshold = (hsize_t)env_size("HDF_ALIGNMENT_THRESHOLD", 1048576); /* 1 MB */
  hsize_t alignment_value     = (hsize_t)env_size("HDF_ALIGNMENT_VALUE", 4194304);     /* 4 MB */

  // hsize_t meta_block_size = (hsize_t)env_size("HDF_META_BLOCK_SIZE", 4194304); /* 4 MB */
  // In THDF_create_zslab_file
  hsize_t meta_block_size = (hsize_t)env_size("HDF_META_BLOCK_SIZE", 4096);

  int coll_metadata_ops   = (int)env_long("HDF_COLL_METADATA_OPS", 1);
  int coll_metadata_write = (int)env_long("HDF_COLL_METADATA_WRITE", 1);

  /* ================================================================
   * Build MPI_Info
   * ================================================================ */
  MPI_Info info;
  MPI_Info_create(&info);

  /* ---- Common hints ---- */
  {
    char buf[64];
    snprintf(buf, sizeof(buf), "%zu", cb_buffer_size);
    MPI_Info_set(info, "cb_buffer_size", buf);
  }
  MPI_Info_set(info, "cb_nodes", cb_nodes_str);
  MPI_Info_set(info, "romio_cb_write", romio_cb_write);
  MPI_Info_set(info, "romio_ds_write", romio_ds_write);
  MPI_Info_set(info, "romio_ds_read", romio_ds_read);

  /* ---- Filesystem-specific hints ---- */
  if (strcmp(fs_type, "lustre") == 0) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%ld", lustre_striping_factor);
    MPI_Info_set(info, "striping_factor", buf);
    snprintf(buf, sizeof(buf), "%zu", lustre_striping_unit);
    MPI_Info_set(info, "striping_unit", buf);
    MPI_Info_set(info, "romio_lustre_start_ioctl_filldir", "0");

  } else if (strcmp(fs_type, "gpfs") == 0) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%ld", gpfs_recv_queues);
    MPI_Info_set(info, "gpfs_recv_queues", buf);
    MPI_Info_set(info, "access_style", gpfs_access_style);
    MPI_Info_set(info, "collective_buffering", "true");
    if (gpfs_num_io_nodes > 0) {
      snprintf(buf, sizeof(buf), "%ld", gpfs_num_io_nodes);
      MPI_Info_set(info, "gpfs_num_io_nodes", buf);
    }

  } else {
    if (mpi_rank == 0)
      fprintf(stderr,
              "[fapl] WARNING: unknown HDF_FS_TYPE='%s', "
              "using common hints only\n",
              fs_type);
  }

  /* ================================================================
   * Build FAPL
   * ================================================================ */
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  if (fapl < 0) {
    MPI_Info_free(&info);
    return -1;
  }

  H5Pset_fapl_mpio(fapl, comm, info);
  H5Pset_alignment(fapl, alignment_threshold, alignment_value);
  H5Pset_meta_block_size(fapl, meta_block_size);

  if (coll_metadata_ops)
    H5Pset_all_coll_metadata_ops(fapl, 1);
  if (coll_metadata_write)
    H5Pset_coll_metadata_write(fapl, 1);

  MPI_Info_free(&info);

  /* ================================================================
   * Print active configuration (rank 0 only)
   * ================================================================ */
  if (mpi_rank == 0) {
    printf("[fapl] %-30s = %s\n", "HDF_FS_TYPE", fs_type);
    printf("[fapl] %-30s = %s\n", "cb_nodes", cb_nodes_str);
    printf("[fapl] %-30s = %zu\n", "cb_buffer_size", cb_buffer_size);
    printf("[fapl] %-30s = %s\n", "romio_cb_write", romio_cb_write);
    printf("[fapl] %-30s = %s\n", "romio_ds_write", romio_ds_write);
    printf("[fapl] %-30s = %zu\n", "alignment_value", (size_t)alignment_value);
    printf("[fapl] %-30s = %zu\n", "meta_block_size", (size_t)meta_block_size);
    printf("[fapl] %-30s = %d\n", "coll_metadata_ops", coll_metadata_ops);
    printf("[fapl] %-30s = %d\n", "coll_metadata_write", coll_metadata_write);

    if (strcmp(fs_type, "lustre") == 0) {
      printf("[fapl] %-30s = %ld\n", "lustre striping_factor", lustre_striping_factor);
      printf("[fapl] %-30s = %zu\n", "lustre striping_unit", lustre_striping_unit);
    } else if (strcmp(fs_type, "gpfs") == 0) {
      printf("[fapl] %-30s = %ld\n", "gpfs recv_queues", gpfs_recv_queues);
      printf("[fapl] %-30s = %s\n", "gpfs access_style", gpfs_access_style);
      if (gpfs_num_io_nodes > 0)
        printf("[fapl] %-30s = %ld\n", "gpfs num_io_nodes", gpfs_num_io_nodes);
    }
  }

  return fapl;
}
