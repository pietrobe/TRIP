
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502984
//           1 23456789 123456789 123456789 1234567
#endif

#include "hdf_output_domain_decomposition_3D.h"
#include "hdf_output_example_3D.h"

void
make_input_example_dataset(double *stokes_IQUI, int N_x, int N_y, int N_z, int N_incl, int N_azimuth, int N_frequencies,
                           int N_stokes, int *stride_x, int *stride_y, int *stride_z, int *stride_incl,
                           int *stride_azimuth, int *stride_frequencies, int *stride_stokes) {
  for (int i = 0; i < N_x; i++) {
    for (int j = 0; j < N_y; j++) {
      for (int k = 0; k < N_z; k++) {
        for (int incl = 0; incl < N_incl; incl++) {
          for (int az = 0; az < N_azimuth; az++) {
            for (int f = 0; f < N_frequencies; f++) {
              for (int s = 0; s < N_stokes; s++) {
                int idx = ((((((i)*N_y + (j)) * N_z + (k)) *  //
                                 N_incl +
                             incl) *
                                N_azimuth +
                            az) *
                               N_frequencies +
                           f) *
                              N_stokes +
                          s;

                double value = 0.0;
                double xx    = (i + 1) / (0.2 * N_x);
                double yy    = (j + 1) / (0.2 * N_y);
                double zz    = (k + 1) / (0.2 * N_z);

                if (s == 0) {  // Stokes I
                  value = sin(xx) * cos(yy + zz) * sin(zz) * (incl + 1) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) *
                          (f + 1) / (0.2 * N_frequencies);
                } else if (s == 1) {  // Stokes Q/I
                  value = sin(2.0 * xx * yy * zz * (incl + 1) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) * (f + 1) /
                              (0.2 * N_frequencies));
                } else if (s == 2) {  // Stokes U/I
                  value = cos(3.0 * xx * yy * zz * (incl + 1) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) * (f + 1) /
                              (0.2 * N_frequencies)) +
                          sin(5.0 * xx * yy * zz * (incl + 1) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) * (f + 1) /
                              (0.2 * N_frequencies));
                } else if (s == 3) {  // Stokes V/I
                  value = cos(4.0 * xx * yy * zz * (incl + 2) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) * (f + 1) /
                              (0.2 * N_frequencies)) -
                          sin(6.0 * xx * yy * zz * (incl + 1) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) * (f + 1) /
                              (0.2 * N_frequencies));
                }

                stokes_IQUI[idx] = value;
              }
            }
          }
        }
      }
    }
  }
  // Calculate strides for the input data layout (C-style contiguous)
  *stride_stokes      = 1;
  *stride_frequencies = N_stokes;
  *stride_azimuth     = N_stokes * N_frequencies;
  *stride_incl        = N_stokes * N_frequencies * N_azimuth;
  *stride_z           = N_stokes * N_frequencies * N_azimuth * N_incl;
  *stride_y           = N_stokes * N_frequencies * N_azimuth * N_incl * N_z;
  *stride_x           = N_stokes * N_frequencies * N_azimuth * N_incl * N_z * N_y;
}

void
make_input_example_dataset_zz(double *stokes_IQUI, int N_x, int N_y, int N_z, int N_incl, int N_azimuth,
                              int N_frequencies, int N_stokes, int *stride_x, int *stride_y, int *stride_z,
                              int *stride_incl, int *stride_azimuth, int *stride_frequencies, int *stride_stokes,
                              int z_par) {
  for (int i = 0; i < N_x; i++) {
    for (int j = 0; j < N_y; j++) {
      for (int k = 0; k < N_z; k++) {
        for (int incl = 0; incl < N_incl; incl++) {
          for (int az = 0; az < N_azimuth; az++) {
            for (int f = 0; f < N_frequencies; f++) {
              for (int s = 0; s < N_stokes; s++) {
                int idx = ((((((i)*N_y + (j)) * N_z + (k)) *  //
                                 N_incl +
                             incl) *
                                N_azimuth +
                            az) *
                               N_frequencies +
                           f) *
                              N_stokes +
                          s;

                double value = 0.0;
                double xx    = (i + 1) / (0.2 * N_x);
                double yy    = (j + 1) / (0.2 * N_y);
                double zz    = (k + 1) / (0.2 * N_z);

                if (s == 0) {  // Stokes I
                  value = z_par + k + sin(xx * f / (0.1 * N_frequencies)) * cos(yy + zz + z_par) +
                          sin(zz + z_par) * (incl + 1) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) * (f + 1) /
                              (0.05 * N_frequencies);
                } else if (s == 1) {  // Stokes Q/I
                  value = sin(2.0 * xx * yy * zz * (incl * (z_par + 1)) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) *
                              (f + 1) / (0.1 * N_frequencies));
                } else if (s == 2) {  // Stokes U/I
                  value = z_par + k +
                          cos(3.0 * xx * yy * zz * (incl + 1) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) * (f + 1) /
                              (0.1 * N_frequencies)) +
                          sin(5.0 * xx * yy * (zz + z_par) * (incl + 1) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) *
                              (f + 1) / (0.1 * N_frequencies));
                } else if (s == 3) {  // Stokes V/I
                  value = z_par + k +
                          cos(4.0 * xx * yy * zz * (incl + 2) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) * (f + 1) /
                              (0.1 * N_frequencies)) -
                          sin(6.0 * xx * yy * zz * (z_par + 1) * (incl + 1) / (0.2 * N_incl) * (az + 1) /
                              (0.2 * N_azimuth) * (f + 1) / (0.1 * N_frequencies));
                }

                stokes_IQUI[idx] = value;
              }
            }
          }
        }
      }
    }
  }
  // Calculate strides for the input data layout (C-style contiguous)
  *stride_stokes      = 1;
  *stride_frequencies = N_stokes;
  *stride_azimuth     = N_stokes * N_frequencies;
  *stride_incl        = N_stokes * N_frequencies * N_azimuth;
  *stride_z           = N_stokes * N_frequencies * N_azimuth * N_incl;
  *stride_y           = N_stokes * N_frequencies * N_azimuth * N_incl * N_z;
  *stride_x           = N_stokes * N_frequencies * N_azimuth * N_incl * N_z * N_y;
}

/////////////////////////////////////////////////////////////////////////
// main_3d_example_v2  —  refactored
//
// Fasi:
//   1. rank 0: crea struttura HDF5 (seriale, nessun FAPL MPI)
//   2. rank 0: scrive metadata (grids, geometry) su file seriale
//              → H5Fclose seriale
//   3. MPI_Barrier — tutti i rank aspettano che rank 0 abbia finito
//   4. tutti: aprono il file con FAPL MPI su persistent_comm
//   5. tutti: THDF_open_3D_field_handler_mpi  (una volta sola)
//   6. loop Z: solo scrittura Stokes, nessun open/close
//   7. tutti: THDF_close_3D_field_handler_mpi + H5Fclose
/////////////////////////////////////////////////////////////////////////
int
main_3d_example_v2(int argc, char **argv) {

  (void)argc;
  (void)argv;

  MPI_Init(&argc, &argv);

  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  /* ----------------------------------------------------------------
   * Parametri
   * ---------------------------------------------------------------- */
  const int  N_x               = read_env_int_or_default("HDF_N_X", 8);
  const int  N_y               = read_env_int_or_default("HDF_N_Y", 8);
  const int  N_z               = read_env_int_or_default("HDF_N_Z", 18);
  const int  N_frequencies     = read_env_int_or_default("HDF_N_FREQUENCIES", 96);
  const int  N_incl            = read_env_int_or_default("HDF_N_INCL", 16);
  const int  N_azimuth         = read_env_int_or_default("HDF_N_AZIMUTH", 8);
  const int  N_stokes          = read_env_int_or_default("HDF_N_STOKES", 4);
  const bool normalized_output = read_env_int_or_default("HDF_NORMALIZE_OUTPUT", 0) != 0;
  int        step_z            = read_env_int_or_default("HDF_STEP_Z", 2);

  /* ----------------------------------------------------------------
   * Allocazione array di griglia
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

  char filename[PATH_MAX];
  build_output_filepath(filename, sizeof(filename), "output_3D_field_mpi.h5");

  /* ================================================================
   * FASE 1+2 — rank 0: crea file + scrive metadata (seriale)
   * ================================================================ */
  if (mpi_rank == 0) {

    const double t0 = MPI_Wtime();

    printf("Rank 0: creating HDF5 structure in '%s'...\n", filename);
    if (THDF_create_3D_field_handler_mpi(filename, normalized_output, N_x, N_y, N_z, N_incl, N_azimuth, N_frequencies) <
        0) {
      fprintf(stderr, "Error creating /radiation_field structure\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    printf("[timing] create structure: %.6f s\n", MPI_Wtime() - t0);

    /* ---- Apri serialmente per scrivere i metadata di griglia ---- */
    hid_t fserial = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    if (fserial < 0) {
      fprintf(stderr, "Rank 0: error opening file for metadata write\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    /* ---- Frequenze ---- */
    THDF_frequencies_grid_t fg;
    fg.N_frequencies = N_frequencies;
    fg.frequencies   = frequencies;
    THDF_write_frequencies_grid_to_hdf5(fserial, &fg);

    /* ---- Griglia angolare ---- */
    THDF_angular_grid_t ag;
    ag.N_azimuthal_angles   = N_azimuth;
    ag.N_inclination_angles = N_incl;
    ag.N_directions         = N_incl * N_azimuth;
    ag.inclination_angles   = theta;
    ag.inclinations_indices = inclinations_idx;
    ag.azimuthal_angles     = chi;
    ag.azimuthal_indices    = azimuthal_idx;
    THDF_write_angular_grid_to_hdf5(fserial, &ag);

    /* ---- Geometria ---- */
    THDF_geometry_3D_t g3;
    g3.N_x     = N_x;
    g3.N_y     = N_y;
    g3.N_z     = N_z;
    g3.heights = heights;
    g3.delta   = 100.0;
    THDF_write_geometry_3D_to_hdf5(fserial, &g3);

    /* Flush esplicito prima del close — garantisce visibilità su GPFS */
    H5Fflush(fserial, H5F_SCOPE_GLOBAL);
    H5Fclose(fserial);

    printf("Rank 0: metadata written, releasing serial file handle\n");
  }

  /* ================================================================
   * FASE 3 — tutti i rank aspettano che rank 0 abbia chiuso il file
   * ================================================================ */
  MPI_Barrier(MPI_COMM_WORLD);

  /* ================================================================
   * Decomposizione dominio
   * ================================================================ */
  int N_local_x, N_local_y, N_local_z;
  int local_start_x, local_start_y, local_start_z;

  decompose_domain_3d(mpi_rank, mpi_size,                              //
                      N_x, N_y, N_z,                                   //
                      &N_local_x, &N_local_y, &N_local_z,              //
                      &local_start_x, &local_start_y, &local_start_z,  //
                      NULL, NULL, NULL);                               //

  if (mpi_rank == 0)
    printf("Rank 0: local domain X[%d+%d] Y[%d+%d] Z[%d+%d]\n", local_start_x, N_local_x, local_start_y, N_local_y,
           local_start_z, N_local_z);

  /* ================================================================
   * FASE 4 — apri il file collettivamente (una volta sola)
   *
   * persistent_comm include SOLO i rank con N_local_z > 0.
   * Con decomposizioni bilanciate coincide con MPI_COMM_WORLD,
   * ma gestisce correttamente il caso edge N_z < mpi_size.
   * ================================================================ */
  const int has_z_data = (N_local_z > 0) ? 1 : MPI_UNDEFINED;

  MPI_Comm persistent_comm = MPI_COMM_NULL;
  MPI_Comm_split(MPI_COMM_WORLD, has_z_data, mpi_rank, &persistent_comm);

  hid_t                    file_id = -1;
  THDF_3D_field_handler_t *handler = NULL;

  if (persistent_comm != MPI_COMM_NULL) {

    /* FAPL con MPI-IO collettivo e hint GPFS */
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "cb_buffer_size", "67108864"); /* 64 MB aggregator buf */
    MPI_Info_set(info, "cb_nodes", "8");

    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, persistent_comm, info);
    H5Pset_alignment(fapl, 1, 4 * 1024 * 1024); /* allinea a 4 MB stripe */
    H5Pset_all_coll_metadata_ops(fapl, 1);
    H5Pset_coll_metadata_write(fapl, 1);
    MPI_Info_free(&info);

    file_id = H5Fopen(filename, H5F_ACC_RDWR, fapl);
    H5Pclose(fapl);

    if (file_id < 0) {
      fprintf(stderr, "Rank %d: error opening file collectively\n", mpi_rank);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    /* ================================================================
     * FASE 5 — apri handler (H5Gopen + H5Dopen × 4-8): UNA VOLTA SOLA
     * ================================================================ */
    handler = THDF_open_3D_field_handler_mpi(file_id, normalized_output, N_x, N_y, N_z, N_incl, N_azimuth, N_frequencies);
    if (!handler) {
      fprintf(stderr, "Rank %d: error opening 3D field handler\n", mpi_rank);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  /* Tutti i rank sincronizzati dopo l'apertura */
  MPI_Barrier(MPI_COMM_WORLD);

  /* ================================================================
   * Allocazione buffer locali
   * ================================================================ */
  int max_N_local_z, min_N_local_z;
  MPI_Allreduce(&N_local_z, &max_N_local_z, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&N_local_z, &min_N_local_z, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  if (max_N_local_z < step_z) {
    if (mpi_rank == 0)
      printf("Warning: step_z (%d) > max_N_local_z (%d), adjusting\n", step_z, max_N_local_z);
    step_z = max_N_local_z;
  }

  if (mpi_rank == 0) {
    printf("Max/min local Z: %d / %d — Z blocks of %d\n", max_N_local_z, min_N_local_z, step_z);
  }

  const size_t input_size  = (size_t)N_local_x * N_local_y * step_z * N_incl * N_azimuth * N_frequencies * N_stokes;
  const size_t output_size = (size_t)N_local_x * N_local_y * step_z * N_incl * N_azimuth * N_frequencies;

  double         *stokes_IQUI = malloc(input_size * sizeof(double));
  THDF_float32_t *stokes_I    = malloc(output_size * sizeof(THDF_float32_t));
  THDF_float32_t *stokes_QI   = malloc(output_size * sizeof(THDF_float32_t));
  THDF_float32_t *stokes_UI   = malloc(output_size * sizeof(THDF_float32_t));
  THDF_float32_t *stokes_VI   = malloc(output_size * sizeof(THDF_float32_t));

  THDF_float_t *norm_I = NULL, *norm_QI = NULL, *norm_UI = NULL, *norm_VI = NULL;
  if (normalized_output) {
    const size_t ns = (size_t)N_local_x * N_local_y * step_z * N_incl * N_azimuth;
    norm_I          = malloc(ns * sizeof(THDF_float_t));
    norm_QI         = malloc(ns * sizeof(THDF_float_t));
    norm_UI         = malloc(ns * sizeof(THDF_float_t));
    norm_VI         = malloc(ns * sizeof(THDF_float_t));
  }

  /* ================================================================
   * FASE 6 — loop Z: solo compute + scrittura, nessun open/close
   * ================================================================ */
  const double write_t0 = MPI_Wtime();

  for (int local_zi = 0; local_zi < max_N_local_z; local_zi += step_z) {

    const int current_nz = (local_zi + step_z <= N_local_z) ? step_z : (N_local_z - local_zi);
    const int is_writer  = (local_zi < N_local_z) ? 1 : 0;
    const int id_z       = local_start_z + local_zi;

    if (is_writer) {

      if (mpi_rank < 4)
        printf("Rank %d: writing Z block [%d:%d]\n", mpi_rank, id_z, id_z + current_nz);

      int sx, sy, sz, si, sa, sf, ss;
      make_input_example_dataset_zz(stokes_IQUI, N_local_x, N_local_y, current_nz, N_incl, N_azimuth, N_frequencies,
                                    N_stokes, &sx, &sy, &sz, &si, &sa, &sf, &ss, id_z);

      THDF_3D_field_t output_field;
      THDF_copy_3D_block_field(&output_field, normalized_output, stokes_IQUI, stokes_I, stokes_QI, stokes_UI, stokes_VI,
                               norm_I, norm_QI, norm_UI, norm_VI, 0, 0, 0, N_local_x, N_local_y, current_nz, N_incl,
                               N_azimuth, N_frequencies, sx, sy, sz, si, sa, sf, ss);

      /* handler già aperto — nessun open/close qui */
      THDF_write_3D_field_dataset_to_hdf5(handler, &output_field, local_start_x, local_start_y, id_z, 0, 0, N_local_x,
                                          N_local_y, current_nz, N_incl, N_azimuth, N_frequencies);
    }

    /*
     * Lockstep tra tutti i rank (inclusi quelli senza dati).
     * Necessario perché persistent_comm potrebbe non includere
     * tutti i rank — la barriera su MPI_COMM_WORLD è l'unico
     * punto di sincronizzazione globale sicuro.
     */
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (mpi_rank == 0)
    printf("[timing] total Z write loop: %.6f s\n", MPI_Wtime() - write_t0);

  /* ================================================================
   * FASE 7 — chiudi handler e file (una volta sola)
   * ================================================================ */
  if (persistent_comm != MPI_COMM_NULL) {
    THDF_close_3D_field_handler_mpi(handler); /* H5Dclose × 4-8, H5Gclose */
    H5Fclose(file_id);                        /* sync GPFS — una volta sola */
    MPI_Comm_free(&persistent_comm);
  }

  /* ================================================================
   * Cleanup
   * ================================================================ */
  free(stokes_IQUI);
  free(stokes_I);
  free(stokes_QI);
  free(stokes_UI);
  free(stokes_VI);
  free(norm_I);
  free(norm_QI);
  free(norm_UI);
  free(norm_VI);
  free(heights);
  free(frequencies);
  free(theta);
  free(chi);
  free(inclinations_idx);
  free(azimuthal_idx);

  MPI_Finalize();
  return 0;
}
