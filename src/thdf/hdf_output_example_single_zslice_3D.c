/**
 * hdf_output_single_file_3D.c
 *
 * Single-file parallel HDF5 writer example using the _single z-slab API.
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
#include "hdf_output_single_zslice_3D.h"
#include "output_utils.h"

/* Set to 0 (e.g. -DZSLAB_PROGRESS_BAR=0) to silence the z-slab progress bar. */
#ifndef ZSLAB_PROGRESS_BAR
#define ZSLAB_PROGRESS_BAR 1
#endif

/* =========================================================================
 * main_3d_example_single_file
 * ========================================================================= */
int
main_3d_example_single_file(int argc, char **argv) {

  (void)argc;
  (void)argv;

  MPI_Init(&argc, &argv);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  if (mpi_rank == 0)
    printf("Running 3D single-file example with MPI: main_3d_example_single_file %s:%d\n", __FILE__, __LINE__);

  const double t_start = MPI_Wtime();

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
  const double total_size = (double)(N_x * N_y * N_z) * (double)(N_incl * N_azimuth * N_frequencies * N_stokes);

  if (mpi_rank == 0) {
    printf("export HDF_N_X=%d\nexport HDF_N_Y=%d\nexport HDF_N_Z=%d\n", N_x, N_y, N_z);
    printf("export HDF_N_FREQUENCIES=%d\nexport HDF_N_INCL=%d\nexport HDF_N_AZIMUTH=%d\n", N_frequencies, N_incl,
           N_azimuth);
    printf("export HDF_N_STOKES=%d\nexport HDF_NORMALIZE_OUTPUT=%d\nexport HDF_STEP_Z=%d\n", N_stokes,
           normalized_output ? 1 : 0, step_z);
    printf("export HDF_OUTPUT_DIR='%s'\n", output_dir);
    printf("[single] set HDF_CHUNK_Z=%d (= HDF_STEP_Z) for zero lock contention\n", step_z);
  }

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

  int N_local_x, N_local_y, N_local_z;
  int local_start_x, local_start_y, local_start_z;
  int Px, Py, Pz;
  decompose_domain_3d(mpi_rank, mpi_size, N_x, N_y, N_z, &N_local_x, &N_local_y, &N_local_z, &local_start_x,
                      &local_start_y, &local_start_z, &Px, &Py, &Pz);
  MPI_Barrier(MPI_COMM_WORLD);

  int max_N_local_z, min_N_local_z;
  MPI_Allreduce(&N_local_z, &max_N_local_z, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&N_local_z, &min_N_local_z, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  char output_file_name[PATH_MAX];
  {
    char fn[PATH_MAX];
    snprintf(fn, sizeof(fn), "radiation_field_3d.h5");
    build_output_filepath(output_file_name, sizeof(output_file_name), fn);
  }

  /* PHASE 1 — create with GLOBAL N_z */
  {
    double t0 = MPI_Wtime();
    if (THDF_create_zslab_file_single(MPI_COMM_WORLD, output_file_name, normalized_output, N_x, N_y, N_z, N_incl,
                                      N_azimuth, N_frequencies) < 0) {
      fprintf(stderr, "Rank %d: create failed\n", mpi_rank);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0)
      printf("[t] create: %.3f s\n", MPI_Wtime() - t0);
  }

  /* PHASE 2 — rank 0 metadata */
  if (mpi_rank == 0) {
    hid_t fs = H5Fopen(output_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
    if (fs < 0) {
      fprintf(stderr, "Rank 0: serial open failed\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    THDF_frequencies_grid_t fg = {N_frequencies, frequencies, HDF_OUT_FLOAT64};
    THDF_write_frequencies_grid_to_hdf5(fs, &fg);
    THDF_angular_grid_t ag = {N_incl * N_azimuth, N_incl,        N_azimuth,      theta, chi,
                              inclinations_idx,   azimuthal_idx, HDF_OUT_FLOAT64};
    THDF_write_angular_grid_to_hdf5(fs, &ag);
    THDF_geometry_3D_t g3 = {N_x, N_y, N_z, heights, 100.0, HDF_OUT_FLOAT64};
    THDF_write_geometry_3D_to_hdf5(fs, &g3);
    H5Fflush(fs, H5F_SCOPE_GLOBAL);
    H5Fclose(fs);
    printf("[t] metadata written\n");
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /* PHASE 3 — open file and handler ONCE with GLOBAL N_z */
  {
    double t0      = MPI_Wtime();
    hid_t  fapl    = build_fapl_from_env(MPI_COMM_WORLD);
    hid_t  file_id = H5Fopen(output_file_name, H5F_ACC_RDWR, fapl);
    H5Pclose(fapl);
    if (file_id < 0) {
      fprintf(stderr, "Rank %d: H5Fopen failed\n", mpi_rank);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    THDF_zslab_handler_t *handler =
        THDF_open_zslab_handler_single(file_id, normalized_output, N_x, N_y, N_z, N_incl, N_azimuth, N_frequencies);
    if (!handler) {
      fprintf(stderr, "Rank %d: open_handler failed\n", mpi_rank);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0)
      printf("[t] open handler: %.3f s\n", MPI_Wtime() - t0);

    /* PHASE 4 — write loop: ALL ranks participate every iteration */
    const size_t    max_out     = (size_t)N_local_x * N_local_y * step_z * N_incl * N_azimuth * N_frequencies;
    const size_t    max_in      = max_out * N_stokes;
    const size_t    max_norm    = (size_t)N_local_x * N_local_y * step_z * N_incl * N_azimuth;
    double         *stokes_IQUI = malloc(max_in * sizeof(double));
    THDF_float32_t *stokes_I    = malloc(max_out * sizeof(THDF_float32_t));
    THDF_float32_t *stokes_QI   = malloc(max_out * sizeof(THDF_float32_t));
    THDF_float32_t *stokes_UI   = malloc(max_out * sizeof(THDF_float32_t));
    THDF_float32_t *stokes_VI   = malloc(max_out * sizeof(THDF_float32_t));
    THDF_float_t   *norm_I      = normalized_output ? malloc(max_norm * sizeof(THDF_float_t)) : NULL;
    THDF_float_t   *norm_QI     = normalized_output ? malloc(max_norm * sizeof(THDF_float_t)) : NULL;
    THDF_float_t   *norm_UI     = normalized_output ? malloc(max_norm * sizeof(THDF_float_t)) : NULL;
    THDF_float_t   *norm_VI     = normalized_output ? malloc(max_norm * sizeof(THDF_float_t)) : NULL;

    const double wt0               = MPI_Wtime();
    const int    total_zslab_steps = (max_N_local_z + step_z - 1) / step_z;

    for (int local_zi = 0; local_zi < max_N_local_z; local_zi += step_z) {
      const int is_writer      = (local_zi < N_local_z);
      const int current_nz     = is_writer ? ((local_zi + step_z <= N_local_z) ? step_z : (N_local_z - local_zi)) : 0;
      const int global_start_z = local_start_z + local_zi;

      THDF_3D_field_t field;
      if (is_writer) {
        int sx, sy, sz, si, sa, sf, ss;
        make_input_example_dataset_zz(stokes_IQUI, N_local_x, N_local_y, current_nz, N_incl, N_azimuth, N_frequencies,
                                      N_stokes, &sx, &sy, &sz, &si, &sa, &sf, &ss, global_start_z);
        THDF_copy_3D_block_field(&field, normalized_output, stokes_IQUI, stokes_I, stokes_QI, stokes_UI, stokes_VI,
                                 norm_I, norm_QI, norm_UI, norm_VI, 0, 0, 0, N_local_x, N_local_y, current_nz, N_incl,
                                 N_azimuth, N_frequencies, sx, sy, sz, si, sa, sf, ss);
      }  // END if is_writer

      double tw0 = MPI_Wtime();

      /* ALL ranks call this — non-writers pass current_nz=0 */
      THDF_write_zslab_block_single(handler, is_writer ? &field : NULL, local_start_x, local_start_y,
                                    global_start_z, /* GLOBAL z offset */
                                    N_local_x, N_local_y, current_nz, N_incl, N_azimuth, N_frequencies);

      MPI_Barrier(MPI_COMM_WORLD);
      if (mpi_rank == 0) {
        printf("[t] z=%d nz=%d: %.3f s\n", global_start_z, current_nz, MPI_Wtime() - tw0);
        const double gb_written = (double)N_x * N_y * current_nz * N_incl * N_azimuth * N_frequencies * N_stokes *
                                  sizeof(THDF_float32_t) / (1024.0 * 1024.0 * 1024.0);
        printf("[t] local written: %.3f GB  approx. peak-bandwidth: %.3f GB/s\n", gb_written,
               gb_written / (MPI_Wtime() - tw0));
      }

#if ZSLAB_PROGRESS_BAR
      const int progess_bar_step = 16;
      if (mpi_rank == 0 && (local_zi / step_z) % progess_bar_step == 0) {
        const int    step_idx = local_zi / step_z + 1;
        const double pct      = (double)step_idx / total_zslab_steps;
        const double elapsed   = MPI_Wtime() - wt0;
        const double eta       = pct > 0.0 ? elapsed / pct - elapsed : 0.0;
        const double gb_so_far = pct * total_size * sizeof(THDF_float32_t) / (1024.0 * 1024.0 * 1024.0);
        const double bw_gbs    = elapsed > 0.0 ? gb_so_far / elapsed : 0.0;
        const int    filled    = (int)(pct * 30);
        char         bar[31];
        for (int b = 0; b < 30; b++) bar[b] = b < filled ? '=' : ' ';
        bar[30] = '\0';
        printf("[prog] [%s] %3.0f%%  step %d/%d  elapsed %.0fs  ETA ~%.0fs  BW ~%.2f GB/s\n", bar, pct * 100.0,
               step_idx, total_zslab_steps, elapsed, eta, bw_gbs);
      }
#endif
    }  // for local_zi

    MPI_Barrier(MPI_COMM_WORLD);

    if (mpi_rank == 0) {
      printf("[t] write loop: %.3f s\n", MPI_Wtime() - wt0);
      double gb_written = (double)total_size * sizeof(THDF_float32_t) / (1024.0 * 1024.0 * 1024.0);
      printf("[t] total written: %.1f GB  effective loop-bandwidth: %.1f GB/s\n", gb_written,
             gb_written / (MPI_Wtime() - wt0));
    }

    free(stokes_IQUI);
    free(stokes_I);
    free(stokes_QI);
    free(stokes_UI);
    free(stokes_VI);
    free(norm_I);
    free(norm_QI);
    free(norm_UI);
    free(norm_VI);

    /* PHASE 5 — close once */
    {
      double t0 = MPI_Wtime();
      THDF_close_zslab_handler_single(handler);
      H5Fflush(file_id, H5F_SCOPE_GLOBAL);
      MPI_Barrier(MPI_COMM_WORLD);
      if (mpi_rank == 0)
        printf("[t] H5Fflush: %.3f s\n", MPI_Wtime() - t0);
      t0 = MPI_Wtime();
      H5Fclose(file_id);
      MPI_Barrier(MPI_COMM_WORLD);
      if (mpi_rank == 0)
        printf("[t] H5Fclose: %.3f s\n", MPI_Wtime() - t0);
    }
  }

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
    printf("[t] total: %.3f s\n", t_total);
    const double gb = total_size * sizeof(THDF_float32_t) / (1024.0 * 1024.0 * 1024.0);
    printf("[sizes] N_x=%d  N_y=%d  N_z=%d  N_incl=%d  N_azimuth=%d  N_frequencies=%d  N_stokes=%d\n", N_x, N_y, N_z,
           N_incl, N_azimuth, N_frequencies, N_stokes);
    printf("[t] total size: %.1f GB\n", gb);
    printf("[MPI] max local z-slab: %d  min local z-slab: %d\n", max_N_local_z, min_N_local_z);
    printf("[MPI] decomposition: (%d, %d, %d) for domain (%d, %d, %d)\n", Px, Py, Pz, N_x, N_y, N_z);
    printf("[MPI] mpi_size: %d\n", mpi_size);
    printf("[t] output: %.1f GB  TTS-bandwidth: %.1f GB/s\n", gb, gb / t_total);
  }

  MPI_Finalize();
  return 0;
}