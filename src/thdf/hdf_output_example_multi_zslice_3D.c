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

#include "hdf_fapl_from_env.h"
#include "hdf_output_domain_decomposition_3D.h"
#include "hdf_output_example_3D.h"
#include "hdf_output_field.h"
#include "hdf_output_zslice_3D.h"
#include "output_utils.h"

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
int
main_3d_example_multi_zslice(int argc, char **argv) {

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
  const int  N_z               = read_env_int_or_default("HDF_N_Z", 13);
  const int  N_frequencies     = read_env_int_or_default("HDF_N_FREQUENCIES", 96);
  const int  N_incl            = read_env_int_or_default("HDF_N_INCL", 16);
  const int  N_azimuth         = read_env_int_or_default("HDF_N_AZIMUTH", 8);
  const int  N_stokes          = read_env_int_or_default("HDF_N_STOKES", 4);
  const bool normalized_output = read_env_int_or_default("HDF_NORMALIZE_OUTPUT", 0) != 0;
  int        step_z            = read_env_int_or_default("HDF_STEP_Z", 2);

  char output_dir[PATH_MAX];
  build_output_filepath(output_dir, sizeof(output_dir), "");

  const double total_size = (double)(N_x * N_y * N_z) *  //
                            (double)(N_incl * N_azimuth * N_frequencies * N_stokes);

  if (mpi_rank == 0) {
    printf("export HDF_N_X=%d\n", N_x);
    printf("export HDF_N_Y=%d\n", N_y);
    printf("export HDF_N_Z=%d\n", N_z);
    printf("export HDF_N_FREQUENCIES=%d\n", N_frequencies);
    printf("export HDF_N_INCL=%d\n", N_incl);
    printf("export HDF_N_AZIMUTH=%d\n", N_azimuth);
    printf("export HDF_N_STOKES=%d\n", N_stokes);
    printf("export HDF_NORMALIZE_OUTPUT=%d\n", normalized_output ? 1 : 0);
    printf("export HDF_STEP_Z=%d\n", step_z);
    printf("export HDF_OUTPUT_DIR='%s'\n", output_dir);
  }

  // /* strip trailing slash from dir */
  // if (output_dir[strlen(output_dir) - 1] == '/')
  //   output_dir[strlen(output_dir) - 1] = '\0';

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

  for (int i = 0; i < N_incl; i++) {
    for (int j = 0; j < N_azimuth; j++) {
      inclinations_idx[i * N_azimuth + j] = i;
      azimuthal_idx[i * N_azimuth + j]    = j;
    }
  }

  /* ----------------------------------------------------------------
   * Domain decomposition
   * ---------------------------------------------------------------- */
  int N_local_x, N_local_y, N_local_z;
  int local_start_x, local_start_y, local_start_z;
  int Px, Py, Pz;

  decompose_domain_3d(mpi_rank, mpi_size,                              //
                      N_x, N_y, N_z,                                   //
                      &N_local_x, &N_local_y, &N_local_z,              //
                      &local_start_x, &local_start_y, &local_start_z,  //
                      &Px, &Py, &Pz);                                  //

  MPI_Barrier(MPI_COMM_WORLD);

  int max_N_local_z, min_N_local_z;
  MPI_Allreduce(&N_local_z, &max_N_local_z, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&N_local_z, &min_N_local_z, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  for (int local_zi = 0; local_zi < max_N_local_z; local_zi += step_z) {

    const int is_writer      = (local_zi < N_local_z);
    const int current_nz     = is_writer ? ((local_zi + step_z <= N_local_z) ? step_z : (N_local_z - local_zi)) : 0;
    const int global_start_z = local_start_z + local_zi;

    /* Colore univoco: evita collisioni tra iterazioni diverse */
    const int split_color = is_writer ? (local_zi * N_z + global_start_z) : MPI_UNDEFINED;

    MPI_Comm slab_comm = MPI_COMM_NULL;
    MPI_Comm_split(MPI_COMM_WORLD, split_color, mpi_rank, &slab_comm);

    if (is_writer) {

      char output_file_name[PATH_MAX];
      {
        char slab_filename[PATH_MAX];
        build_slab_filename_slice(slab_filename, sizeof(slab_filename), global_start_z, current_nz, 0, N_x, 0, N_y);
        build_output_filepath(output_file_name, sizeof(output_file_name), slab_filename);
      }

      int slab_main_rank;
      MPI_Allreduce(&mpi_rank, &slab_main_rank, 1, MPI_INT, MPI_MIN, slab_comm);

      // printf("Slab z: %d main rank %d\n", global_start_z, slab_main_rank);

      const int err_create = THDF_create_zslab_file(slab_comm, output_file_name,        //
                                                    normalized_output,                  //
                                                    N_x, N_y, current_nz,               //
                                                    N_incl, N_azimuth, N_frequencies);  //
      if (err_create < 0) {
        fprintf(stderr, "Rank %d: failed to create slab file\n", mpi_rank);
        MPI_Abort(MPI_COMM_WORLD, -1);
      }

      MPI_Barrier(slab_comm);

      if (mpi_rank == slab_main_rank) {
        hid_t fs = H5Fopen(output_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
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
        THDF_geometry_3D_t g3 = {N_x, N_y, current_nz, heights + global_start_z, 100.0};
        THDF_write_geometry_3D_to_hdf5(fs, &g3);

        H5Fflush(fs, H5F_SCOPE_GLOBAL);
        H5Fclose(fs);
      }

      MPI_Barrier(slab_comm);

      {  // Last phase: wrting data
         // MPI_Info info;
         // MPI_Info_create(&info);
         // MPI_Info_set(info, "cb_buffer_size", "134217728"); /* 128 MB */
         // MPI_Info_set(info, "cb_nodes", "8");

        // hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        // H5Pset_fapl_mpio(fapl, slab_comm, info);
        // H5Pset_alignment(fapl, 1, 4 * 1024 * 1024);
        // H5Pset_all_coll_metadata_ops(fapl, 1);
        // H5Pset_coll_metadata_write(fapl, 1);
        // MPI_Info_free(&info);

        hid_t fapl = build_fapl_from_env(slab_comm);

        hid_t file_id = H5Fopen(output_file_name, H5F_ACC_RDWR, fapl);
        H5Pclose(fapl);
        if (file_id < 0) {
          fprintf(stderr, "Rank %d: collective H5Fopen failed\n", mpi_rank);
          MPI_Abort(MPI_COMM_WORLD, -1);
        }

        double t0 = MPI_Wtime();

        THDF_zslab_handler_t *handler =
            THDF_open_zslab_handler(file_id, normalized_output, N_x, N_y, current_nz,  // ← fix principale
                                    N_incl, N_azimuth, N_frequencies);
        if (!handler) {
          fprintf(stderr, "Rank %d: THDF_open_zslab_handler failed\n", mpi_rank);
          MPI_Abort(MPI_COMM_WORLD, -1);
        }

        MPI_Barrier(slab_comm);
        if (mpi_rank == 0)
          printf("[t] slab %05d open handler: %.3f s\n", global_start_z, MPI_Wtime() - t0);

        const size_t output_size = (size_t)N_local_x * N_local_y * current_nz * N_incl * N_azimuth * N_frequencies;
        const size_t input_size  = output_size * N_stokes;
        const size_t norm_size   = (size_t)N_local_x * N_local_y * current_nz * N_incl * N_azimuth;

        double         *stokes_IQUI = malloc(input_size * sizeof(double));
        THDF_float32_t *stokes_I    = malloc(output_size * sizeof(THDF_float32_t));
        THDF_float32_t *stokes_QI   = malloc(output_size * sizeof(THDF_float32_t));
        THDF_float32_t *stokes_UI   = malloc(output_size * sizeof(THDF_float32_t));
        THDF_float32_t *stokes_VI   = malloc(output_size * sizeof(THDF_float32_t));
        THDF_float_t   *norm_I      = normalized_output ? malloc(norm_size * sizeof(THDF_float_t)) : NULL;
        THDF_float_t   *norm_QI     = normalized_output ? malloc(norm_size * sizeof(THDF_float_t)) : NULL;
        THDF_float_t   *norm_UI     = normalized_output ? malloc(norm_size * sizeof(THDF_float_t)) : NULL;
        THDF_float_t   *norm_VI     = normalized_output ? malloc(norm_size * sizeof(THDF_float_t)) : NULL;

        int sx, sy, sz, si, sa, sf, ss;
        make_input_example_dataset_zz(stokes_IQUI, N_local_x, N_local_y, current_nz, N_incl, N_azimuth, N_frequencies,
                                      N_stokes, &sx, &sy, &sz, &si, &sa, &sf, &ss, global_start_z);

        THDF_3D_field_t field;
        THDF_copy_3D_block_field(&field, normalized_output, stokes_IQUI, stokes_I, stokes_QI, stokes_UI, stokes_VI,
                                 norm_I, norm_QI, norm_UI, norm_VI, 0, 0, 0, N_local_x, N_local_y, current_nz, N_incl,
                                 N_azimuth, N_frequencies, sx, sy, sz, si, sa, sf, ss);

        const double tw0 = MPI_Wtime();
        THDF_write_zslab_block(handler, &field,               //
                               local_start_x, local_start_y,  //
                               0,                             /* z offset in slab file — 0-based */
                               N_local_x, N_local_y, current_nz, N_incl, N_azimuth, N_frequencies);
        const double tw0_clock = MPI_Wtime() - tw0;
        if (mpi_rank == 0)
          printf("[t] write time: %f\n", tw0_clock);

        MPI_Barrier(slab_comm);

        ////////////////////////////////////
        THDF_close_zslab_handler(handler);
        H5Fclose(file_id);

        // MPI_Barrier(slab_comm);
        // MPI_Comm_free(&slab_comm);

        free(stokes_IQUI);
        free(stokes_I);
        free(stokes_QI);
        free(stokes_UI);
        free(stokes_VI);
        free(norm_I);
        free(norm_QI);
        free(norm_UI);
        free(norm_VI);
      }  // END IF WRITER
    }

    /* All ranks must hit the same world barrier sequence, even when a rank has
     * no work for this z-step. */
    MPI_Barrier(MPI_COMM_WORLD);

    if (slab_comm != MPI_COMM_NULL)
      MPI_Comm_free(&slab_comm);

    MPI_Barrier(MPI_COMM_WORLD);
  }

  //// free resources
  free(heights);
  free(frequencies);
  free(theta);
  free(chi);
  free(inclinations_idx);
  free(azimuthal_idx);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  const double t_total = MPI_Wtime() - t_start;
  if (mpi_rank == 0) {
    printf("[t] total execution time: %.3f s\n", t_total);
    double total_gb = total_size * sizeof(THDF_float32_t) / (1024 * 1024 * 1024.0);
    printf("[t] total output size: %.1f GB\n", total_gb);
    printf("[t] effective write bandwidth: %.1f GB/s\n", total_gb / t_total);
  }

  MPI_Finalize();

  return 0;
}