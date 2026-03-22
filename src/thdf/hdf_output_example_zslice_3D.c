/**
 * hdf5_zslab_writer.c
 *
 * Parallel HDF5 writer — z-slab strategy.
 *
 * Layout:
 *   slab_z00000.h5  →  [N_x, N_y, N_local_z, N_incl, N_azimuth, N_freq]
 *   slab_z00001.h5  →  [N_x, N_y, N_local_z, ...]
 *   ...
 *   master.h5       →  VDS  [N_x, N_y, N_z, N_incl, N_azimuth, N_freq]
 *
 * Each slab file owns the FULL xy plane for its z range.
 * All ranks sharing the same local_start_z write collectively to one slab.
 * Slab files are written INDEPENDENTLY — zero inter-slab coordination.
 *
 * Design decisions:
 *   - ALLOC_TIME_EARLY removed  (was pre-allocating 385 GB from rank 0)
 *   - Chunk size ~12 MB         (was 768 MB)
 *   - Handler opened once       (was opened per Z block)
 *   - Collective DXPL           (H5FD_MPIO_COLLECTIVE in write_stokes_hyperslab)
 */

#include <errno.h>
#include <hdf5.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "hdf_output_domain_decomposition_3D.h"
#include "hdf_output_example_3D.h"
#include "hdf_output_zslice_3D.h"
#include "output_utils.h"

/* =========================================================================
 * main_3d_example_v2  —  z-slab strategy
 * =========================================================================
 *
 * Phase overview:
 *
 *  1. All ranks  → THDF_create_zslab_file()   (collective on slab_comm)
 *                  34 slab files created in PARALLEL, no inter-slab sync
 *
 *  2. slab_rank 0 → H5Fopen (serial) → write freq/angular/geometry → H5Fclose
 *     MPI_Barrier(slab_comm)
 *
 *  3. All ranks  → H5Fopen (collective FAPL)
 *                  THDF_open_zslab_handler()   ONE TIME
 *
 *  4. Z block loop:
 *       compute data
 *       THDF_write_zslab_block()   (collective DXPL, no open/close)
 *       MPI_Barrier(MPI_COMM_WORLD)
 *
 *  5. THDF_close_zslab_handler() + H5Fclose()   ONE TIME
 *
 *  6. rank 0 → create_zslab_vds_master()
 */
int
main_3d_example_zslice(int argc, char **argv) {

  (void)argc;
  (void)argv;

  MPI_Init(&argc, &argv);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  const double t_start = MPI_Wtime();

  /* ----------------------------------------------------------------
   * Parameters
   * ---------------------------------------------------------------- */
  const int  N_x               = read_env_int_or_default("HDF_N_X", 8);
  const int  N_y               = read_env_int_or_default("HDF_N_Y", 8);
  const int  N_z               = read_env_int_or_default("HDF_N_Z", 12);
  const int  N_frequencies     = read_env_int_or_default("HDF_N_FREQUENCIES", 96);
  const int  N_incl            = read_env_int_or_default("HDF_N_INCL", 16);
  const int  N_azimuth         = read_env_int_or_default("HDF_N_AZIMUTH", 8);
  const int  N_stokes          = read_env_int_or_default("HDF_N_STOKES", 4);
  const bool normalized_output = read_env_int_or_default("HDF_NORMALIZE_OUTPUT", 0) != 0;
  int        step_z            = read_env_int_or_default("HDF_STEP_Z", 1);

  char output_dir[PATH_MAX];
  build_output_filepath(output_dir, sizeof(output_dir), "");
  /* strip trailing slash from dir */

  /* ----------------------------------------------------------------
   * Grid arrays
   * ---------------------------------------------------------------- */
  double *frequencies      = malloc(N_frequencies * sizeof(double));
  double *theta            = malloc(N_incl * sizeof(double));
  double *chi              = malloc(N_azimuth * sizeof(double));
  int    *inclinations_idx = malloc(N_incl * N_azimuth * sizeof(int));
  int    *azimuthal_idx    = malloc(N_incl * N_azimuth * sizeof(int));
  double *heights          = malloc(N_z * sizeof(double));

  for (int i = 0; i < N_frequencies; i++) frequencies[i] = 1.0e14 + i * 1.0e12;
  for (int k = 0; k < N_z; k++) heights[k] = k * 100.0;
  for (int i = 0; i < N_incl; i++) theta[i] = (M_PI / (N_incl - 1)) * i;
  for (int j = 0; j < N_azimuth; j++) chi[j] = j * (2.0 * M_PI / N_azimuth);
  for (int i = 0; i < N_incl; i++)
    for (int j = 0; j < N_azimuth; j++) {
      inclinations_idx[i * N_azimuth + j] = i;
      azimuthal_idx[i * N_azimuth + j]    = j;
    }

  /* ----------------------------------------------------------------
   * Domain decomposition
   * ---------------------------------------------------------------- */
  int N_local_x, N_local_y, N_local_z;
  int local_start_x, local_start_y, local_start_z;

  decompose_domain_3d(mpi_rank, mpi_size,                              //
                      N_x, N_y, N_z,                                   //
                      &N_local_x, &N_local_y, &N_local_z,              //
                      &local_start_x, &local_start_y, &local_start_z,  //
                      NULL, NULL, NULL);                               //

  /* ================================================================
   * Z-slab topology
   * ================================================================ */
  THDF_zslab_topology_t topo = build_zslab_topology(MPI_COMM_WORLD, local_start_z, N_local_z, mpi_rank);

  if (mpi_rank == 0)
    printf("Z-slab topology: %d slabs × %d ranks/slab\n", topo.n_slabs, topo.slab_size);

  char slab_filename[PATH_MAX];
  build_zslab_filename(slab_filename, sizeof(slab_filename), output_dir, topo.slab_id);

  /*
   * Slab file size per dataset:
   *   N_x × N_y × N_local_z × N_incl × N_azimuth × N_freq × 4 byte
   *   = 128 × 128 × (128/n_slabs) × 16 × 8 × 96 × 4
   * With 34 slabs: ~3 GB per file  ← safe for ALLOC_TIME_EARLY if desired
   */

  /* ================================================================
   * PHASE 1 — create slab file structure (collective on slab_comm)
   *
   * All 34 slabs create their file SIMULTANEOUSLY.
   * No inter-slab coordination — pure parallel I/O on GPFS.
   * ================================================================ */
  {
    double t0 = MPI_Wtime();
    if (THDF_create_zslab_file(topo.slab_comm, slab_filename, normalized_output, N_x, N_y, N_local_z, N_incl, N_azimuth,
                               N_frequencies) < 0) {
      fprintf(stderr, "Rank %d: failed to create slab file\n", mpi_rank);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    MPI_Barrier(topo.slab_comm);
    if (topo.slab_rank == 0)
      printf("[t] slab %05d create: %.3f s\n", topo.slab_id, MPI_Wtime() - t0);
  }

  /* ================================================================
   * PHASE 2 — slab_rank 0: write grid metadata (serial)
   *
   * Each slab writes its own metadata independently.
   * heights[] is sliced to the local z range.
   * ================================================================ */
  if (topo.slab_rank == 0) {
    hid_t fs = H5Fopen(slab_filename, H5F_ACC_RDWR, H5P_DEFAULT);
    if (fs < 0) {
      fprintf(stderr, "Rank %d: serial open failed for metadata\n", mpi_rank);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    /* Frequencies and angles: same for every slab */
    THDF_frequencies_grid_t fg = {N_frequencies, frequencies};
    THDF_write_frequencies_grid_to_hdf5(fs, &fg);

    THDF_angular_grid_t ag = {N_incl * N_azimuth, N_incl, N_azimuth, theta, chi, inclinations_idx, azimuthal_idx};
    THDF_write_angular_grid_to_hdf5(fs, &ag);

    /*
     * Geometry: N_z = N_local_z (slab thickness).
     * heights[] points into the global array at the slab's z start.
     */
    THDF_geometry_3D_t g3 = {N_x, N_y, N_local_z, heights + topo.slab_id, 100.0};
    THDF_write_geometry_3D_to_hdf5(fs, &g3);

    H5Fflush(fs, H5F_SCOPE_GLOBAL);
    H5Fclose(fs);
  }
  /* All ranks in the slab wait for rank 0 to finish */
  MPI_Barrier(topo.slab_comm);

  /* ================================================================
   * PHASE 3 — open persistent handler (once for the entire job)
   * ================================================================ */
  MPI_Info info;
  MPI_Info_create(&info);
  MPI_Info_set(info, "cb_buffer_size", "134217728"); /* 128 MB */
  MPI_Info_set(info, "cb_nodes", "8");

  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl, topo.slab_comm, info);
  H5Pset_alignment(fapl, 1, 4 * 1024 * 1024);
  H5Pset_all_coll_metadata_ops(fapl, 1);
  H5Pset_coll_metadata_write(fapl, 1);
  MPI_Info_free(&info);

  hid_t file_id = H5Fopen(slab_filename, H5F_ACC_RDWR, fapl);
  H5Pclose(fapl);
  if (file_id < 0) {
    fprintf(stderr, "Rank %d: collective H5Fopen failed\n", mpi_rank);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  {
    double t0 = MPI_Wtime();

    THDF_zslab_handler_t *handler = THDF_open_zslab_handler(file_id, normalized_output,         //
                                                            N_x, N_y, N_local_z,                //
                                                            N_incl, N_azimuth, N_frequencies);  //
    if (!handler) {
      fprintf(stderr, "Rank %d: THDF_open_zslab_handler failed\n", mpi_rank);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    MPI_Barrier(topo.slab_comm);
    if (topo.slab_rank == 0)
      printf("[t] slab %05d open handler: %.3f s\n", topo.slab_id, MPI_Wtime() - t0);

    /* ================================================================
     * PHASE 4 — Z block loop
     *
     * Each iteration writes step_z z-planes to the slab file.
     *
     * Offset rules:
     *   x in file = local_start_x   (global — file spans full N_x)
     *   y in file = local_start_y   (global — file spans full N_y)
     *   z in file = local_zi        (0-based in slab — NOT global z)
     * ================================================================ */
    int max_N_local_z, min_N_local_z;
    MPI_Allreduce(&N_local_z, &max_N_local_z, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&N_local_z, &min_N_local_z, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if (step_z <= 0) {
      if (mpi_rank == 0)
        fprintf(stderr, "Invalid HDF_STEP_Z=%d; forcing step_z=1 to avoid infinite loop\n", step_z);
      step_z = 1;
    }

    if (step_z > max_N_local_z)
      step_z = max_N_local_z;

    if (mpi_rank == 0)
      printf("Z-block loop: step_z=%d  max_local_z=%d  min_local_z=%d\n", step_z, max_N_local_z, min_N_local_z);

    const size_t output_size = (size_t)N_local_x * N_local_y * step_z * N_incl * N_azimuth * N_frequencies;
    const size_t input_size  = output_size * N_stokes;
    const size_t norm_size   = (size_t)N_local_x * N_local_y * step_z * N_incl * N_azimuth;

    double         *stokes_IQUI = malloc(input_size * sizeof(double));
    THDF_float32_t *stokes_I    = malloc(output_size * sizeof(THDF_float32_t));
    THDF_float32_t *stokes_QI   = malloc(output_size * sizeof(THDF_float32_t));
    THDF_float32_t *stokes_UI   = malloc(output_size * sizeof(THDF_float32_t));
    THDF_float32_t *stokes_VI   = malloc(output_size * sizeof(THDF_float32_t));
    THDF_float_t   *norm_I      = normalized_output ? malloc(norm_size * sizeof(THDF_float_t)) : NULL;
    THDF_float_t   *norm_QI     = normalized_output ? malloc(norm_size * sizeof(THDF_float_t)) : NULL;
    THDF_float_t   *norm_UI     = normalized_output ? malloc(norm_size * sizeof(THDF_float_t)) : NULL;
    THDF_float_t   *norm_VI     = normalized_output ? malloc(norm_size * sizeof(THDF_float_t)) : NULL;

    const double write_t0 = MPI_Wtime();

    for (int local_zi = 0; local_zi < max_N_local_z; local_zi += step_z) {

      const int current_nz = (local_zi + step_z <= N_local_z) ? step_z : (N_local_z - local_zi);
      const int is_writer  = (local_zi < N_local_z);
      const int write_nz   = is_writer ? current_nz : 0;

      /*
       * id_z_global: used by make_input to generate physically correct data
       *              (e.g. height-dependent source functions)
       * local_zi:    z offset INSIDE the slab file (always 0-based)
       */
      const int id_z_global = local_start_z + local_zi;

      double tw0 = 0.0;

      if (is_writer) {

        if (mpi_rank < 4)
          printf(
              "Rank %d: writing z block [%d:%d] "
              "(file offset z=%d)\n",
              mpi_rank, id_z_global, id_z_global + current_nz, local_zi);

        int sx, sy, sz, si, sa, sf, ss;
        make_input_example_dataset_zz(stokes_IQUI, N_local_x, N_local_y, current_nz, N_incl, N_azimuth, N_frequencies,
                                      N_stokes, &sx, &sy, &sz, &si, &sa, &sf, &ss, id_z_global);

        THDF_3D_field_t field;
        THDF_copy_3D_block_field(&field, normalized_output, stokes_IQUI, stokes_I, stokes_QI, stokes_UI, stokes_VI,
                                 norm_I, norm_QI, norm_UI, norm_VI, 0, 0, 0, N_local_x, N_local_y, current_nz, N_incl,
                                 N_azimuth, N_frequencies, sx, sy, sz, si, sa, sf, ss);

        tw0 = MPI_Wtime();
        THDF_write_zslab_block(handler, &field, local_start_x, local_start_y,
                               local_zi, /* z offset in slab file — 0-based */
                               N_local_x, N_local_y, write_nz, N_incl, N_azimuth, N_frequencies);
      } else {
        THDF_3D_field_t empty_field = {0};
        THDF_write_zslab_block(handler, &empty_field, local_start_x, local_start_y,
                               local_zi, /* z offset in slab file — 0-based */
                               N_local_x, N_local_y, write_nz, N_incl, N_azimuth, N_frequencies);
      }

      MPI_Barrier(topo.slab_comm);

      if (is_writer && topo.slab_rank == 0) {
        double dt  = MPI_Wtime() - tw0;
        double gbs = (double)N_x * N_y * current_nz * N_incl * N_azimuth * N_frequencies * 4 * 4 / 1e9 /
                     dt; /* 4 datasets × 4 byte */
        printf("[t] slab %05d z=%d write: %.3f s  %.1f GB/s\n", topo.slab_id, local_zi, dt, gbs);
      }

      /*
       * Global barrier: keeps all slabs in lockstep.
       * Not strictly required for correctness (slabs are independent)
       * but prevents fast slabs from overloading GPFS MDT while
       * slow slabs are still writing.
       */
      MPI_Barrier(MPI_COMM_WORLD);
    }

    if (mpi_rank == 0)
      printf("[t] total write loop: %.3f s\n", MPI_Wtime() - write_t0);

    free(stokes_IQUI);
    free(stokes_I);
    free(stokes_QI);
    free(stokes_UI);
    free(stokes_VI);
    free(norm_I);
    free(norm_QI);
    free(norm_UI);
    free(norm_VI);

    /* ================================================================
     * PHASE 5 — close handler and file (once)
     * ================================================================ */
    THDF_close_zslab_handler(handler); /* H5Dclose×8, H5Gclose */
  }
  H5Fclose(file_id); /* collective flush to GPFS */

  MPI_Comm_free(&topo.slab_comm);

  /* ================================================================
   * PHASE 6 — rank 0 builds the VDS master
   *
   * Gathers slab_z_start and N_local_z from every slab root,
   * then writes master.h5 pointing to all slab files.
   * ================================================================ */
  MPI_Barrier(MPI_COMM_WORLD);

  if (topo.slab_rank == 0) {
    int *all_z_start = NULL;
    int *all_nz      = NULL;

    if (mpi_rank == 0) {
      all_z_start = malloc(topo.n_slabs * sizeof(int));
      all_nz      = malloc(topo.n_slabs * sizeof(int));
    }

    MPI_Gather(&topo.slab_id, 1, MPI_INT, all_z_start, 1, MPI_INT, 0, topo.root_comm);
    MPI_Gather(&N_local_z, 1, MPI_INT, all_nz, 1, MPI_INT, 0, topo.root_comm);

    if (mpi_rank == 0) {
      double t0 = MPI_Wtime();
      char   master_path[PATH_MAX];
      snprintf(master_path, sizeof(master_path), "%s", output_dir);

      create_zslab_vds_master(master_path, "master.h5", topo.n_slabs, N_x, N_y, N_z, N_incl, N_azimuth, N_frequencies,
                              normalized_output, all_z_start, all_nz);
      printf("[t] VDS master: %.3f s  →  %s/master.h5\n", MPI_Wtime() - t0, master_path);
      free(all_z_start);
      free(all_nz);
    }

    MPI_Comm_free(&topo.root_comm);
  }

  /* ----------------------------------------------------------------
   * Cleanup
   * ---------------------------------------------------------------- */
  free(heights);
  free(frequencies);
  free(theta);
  free(chi);
  free(inclinations_idx);
  free(azimuthal_idx);

  MPI_Barrier(MPI_COMM_WORLD);
  const double t_total = MPI_Wtime() - t_start;
  if (mpi_rank == 0) {
    printf("[t] total execution time: %.3f s\n", t_total);
    double total_gb = (double)(N_x * N_y * N_z) * (double)(N_incl * N_azimuth * N_frequencies) * 4 *
                      sizeof(THDF_float_t) / (1024 * 1024 * 1024.0);
    printf("[t] total output size: %.1f GB\n", total_gb);
    printf("[t] effective write bandwidth: %.1f GB/s\n", total_gb / t_total);
  }

  MPI_Finalize();
  return 0;
}
