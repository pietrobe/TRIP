
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf_output_example_3D.h"

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

static int
read_env_int_or_default(const char *name, int default_value) {
  const char *value = getenv(name);
  if (value == NULL || *value == '\0') {
    return default_value;
  }

  char *end    = NULL;
  long  parsed = strtol(value, &end, 10);
  if (end == value || *end != '\0' || parsed <= 0 || parsed > INT_MAX) {
    return default_value;
  }

  return (int)parsed;
}

static void
build_output_filepath(char *out_path, size_t out_path_size, const char *default_filename) {
  const char *base_dir = getenv("HDF_OUT_PATH");
  if (base_dir == NULL || *base_dir == '\0') {
    (void)snprintf(out_path, out_path_size, "%s", default_filename);
    return;
  }

  const size_t base_len = strlen(base_dir);
  const int    has_sep  = (base_len > 0 && base_dir[base_len - 1] == '/');

  int written = snprintf(out_path, out_path_size, has_sep ? "%s%s" : "%s/%s", base_dir, default_filename);
  if (written < 0 || (size_t)written >= out_path_size) {
    fprintf(stderr, "[hdf_output_example_3D] WARNING: output path too long, using '%s'\n", default_filename);
    (void)snprintf(out_path, out_path_size, "%s", default_filename);
  }
}

/**
 * Factorizes 'size' into (Px, Py, Pz) such that Px*Py*Pz == size
 * and the processor grid matches the domain aspect ratio as closely as possible.
 */
void
factorize_procs(int size, int Nx, int Ny, int Nz, int *Px, int *Py, int *Pz) {
  int    best_px = 1, best_py = 1, best_pz = size;
  double best_score = 1e18;

  for (int px = 1; px <= size; px++) {
    if (size % px != 0)
      continue;
    for (int py = 1; py <= size / px; py++) {
      if ((size / px) % py != 0)
        continue;
      int pz = size / (px * py);

      /* Local subdomain extents for this factorization */
      double lx = (double)Nx / px;
      double ly = (double)Ny / py;
      double lz = (double)Nz / pz;

      /* Score: sum of pairwise aspect ratios (minimum = 6.0 for a cube) */
      double score = (lx / ly + ly / lx) + (ly / lz + lz / ly) + (lx / lz + lz / lx);

      // printf("  (%d x %d x %d)  local=(%.1f x %.1f x %.1f)  score=%.3f%s\n", px, py, pz, lx, ly, lz, score,
      //        score < best_score ? "  <-- new best" : "");

      if (score < best_score) {
        best_score = score;
        best_px    = px;
        best_py    = py;
        best_pz    = pz;
      }
    }
  }

  *Px = best_px;
  *Py = best_py;
  *Pz = best_pz;
}

/**
 * Distributes n_global points among p_dims ranks.
 * First `remainder` ranks get (base+1), the rest get base.
 */
static void
get_1d_decomposition(int n_global, int p_dims, int p_coord, int *n_local, int *start) {
  int base      = n_global / p_dims;
  int remainder = n_global % p_dims;

  if (p_coord < remainder) {
    *n_local = base + 1;
    *start   = p_coord * (base + 1);
  } else {
    *n_local = base;
    *start   = remainder * (base + 1) + (p_coord - remainder) * base;
  }
}

/**
 * 3D Domain Decomposition — no MPI Cartesian communicator.
 *
 * Maps each MPI rank to a unique (ix, iy, iz) coordinate in a
 * Px x Py x Pz processor grid using simple row-major indexing:
 *
 *   rank = ix * (Py * Pz) + iy * Pz + iz
 *
 * Inputs:
 *   N_x, N_y, N_z  — global domain extents
 *   comm            — MPI communicator
 *
 * Outputs:
 *   N_local_{x,y,z}      — size of this rank's subdomain
 *   local_start_{x,y,z}  — global start index of this rank's subdomain
 *   P{x,y,z}             — processor grid dimensions (optional, can be NULL)
 */
void
decompose_domain_3d(int N_x, int N_y, int N_z, MPI_Comm comm, int *N_local_x, int *N_local_y, int *N_local_z,
                    int *local_start_x, int *local_start_y, int *local_start_z, int *out_Px, int *out_Py, int *out_Pz) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  /* --- 1. Find best (Px, Py, Pz) factorization of 'size' --- */
  int Px, Py, Pz;
  factorize_procs(size, N_x, N_y, N_z, &Px, &Py, &Pz);

  if (rank == 0)
    printf("[decompose_domain_3d] grid: %d x %d x %d procs over %d x %d x %d domain\n", Px, Py, Pz, N_x, N_y, N_z);

  /* --- 2. Map linear rank -> 3D processor coordinate (row-major) ---
   *
   *   rank = ix*(Py*Pz) + iy*Pz + iz
   */
  int ix = rank / (Py * Pz);
  int iy = (rank % (Py * Pz)) / Pz;
  int iz = rank % Pz;

  /* --- 3. Compute local subdomain size and global start offset --- */
  get_1d_decomposition(N_x, Px, ix, N_local_x, local_start_x);
  get_1d_decomposition(N_y, Py, iy, N_local_y, local_start_y);
  get_1d_decomposition(N_z, Pz, iz, N_local_z, local_start_z);

  /* --- 4. Expose grid dims if caller needs them (e.g. to find neighbors) --- */
  if (out_Px)
    *out_Px = Px;
  if (out_Py)
    *out_Py = Py;
  if (out_Pz)
    *out_Pz = Pz;
}

/**
 * Given the processor grid (Px, Py, Pz), returns the rank of a neighbor.
 * Returns MPI_PROC_NULL for out-of-bounds neighbors (domain boundary).
 *
 * Usage example — find the +X neighbor:
 *   int nbr = get_neighbor_rank(Px, Py, Pz, ix+1, iy, iz);
 */
int
get_neighbor_rank(int Px, int Py, int Pz, int ix, int iy, int iz) {
  if (ix < 0 || ix >= Px)
    return MPI_PROC_NULL;
  if (iy < 0 || iy >= Py)
    return MPI_PROC_NULL;
  if (iz < 0 || iz >= Pz)
    return MPI_PROC_NULL;
  return ix * (Py * Pz) + iy * Pz + iz;
}

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
                  value = sin(xx * f / (0.1 * N_frequencies)) * cos(yy + zz + z_par) +
                          sin(zz + z_par) * (incl + 1) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) * (f + 1) /
                              (0.05 * N_frequencies);
                } else if (s == 1) {  // Stokes Q/I
                  value = sin(2.0 * xx * yy * zz * (incl * (z_par + 1)) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) *
                              (f + 1) / (0.1 * N_frequencies));
                } else if (s == 2) {  // Stokes U/I
                  value = cos(3.0 * xx * yy * zz * (incl + 1) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) * (f + 1) /
                              (0.1 * N_frequencies)) +
                          sin(5.0 * xx * yy * (zz + z_par) * (incl + 1) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) *
                              (f + 1) / (0.1 * N_frequencies));
                } else if (s == 3) {  // Stokes V/I
                  value = cos(4.0 * xx * yy * zz * (incl + 2) / (0.2 * N_incl) * (az + 1) / (0.2 * N_azimuth) * (f + 1) /
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

int
main_3d_example(int argc, char **argv) {

  (void)argc;
  (void)argv;

  MPI_Init(&argc, &argv);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int  N_x               = read_env_int_or_default("HDF_N_X", 32);
  int  N_y               = read_env_int_or_default("HDF_N_Y", 32);
  int  N_z               = read_env_int_or_default("HDF_N_Z", 4);
  int  N_frequencies     = read_env_int_or_default("HDF_N_FREQUENCIES", 96);
  int  N_incl            = read_env_int_or_default("HDF_N_INCL", 8);
  int  N_azimuth         = read_env_int_or_default("HDF_N_AZIMUTH", 16);
  int  N_stokes          = read_env_int_or_default("HDF_N_STOKES", 4);
  bool normalized_output = read_env_int_or_default("HDF_NORMALIZE_OUTPUT", 0) != 0;

  double *frequencies = (double *)malloc(N_frequencies * sizeof(double));

  for (int i = 0; i < N_frequencies; i++) {
    frequencies[i] = 1.0e14 + i * 1.0e12;
  }

  double *theta                = (double *)malloc(N_incl * sizeof(double));
  int    *inclinations_indices = (int *)malloc(N_incl * sizeof(int));

  double *chi               = (double *)malloc(N_azimuth * sizeof(double));
  int    *azimuthal_indices = (int *)malloc(N_azimuth * sizeof(int));

  for (int i = 0; i < N_incl; i++) {
    theta[i]                = (i) * (3.141592653589793 / N_incl);  // just example values in radians
    inclinations_indices[i] = i;
  }

  for (int j = 0; j < N_azimuth; j++) {
    chi[j]               = j * (3.141592653589793 * 2.0 / N_azimuth);  // just example values
    azimuthal_indices[j] = j;
  }

  char filename[PATH_MAX];
  build_output_filepath(filename, sizeof(filename), "output_3D_field_mpi.h5");

  if (mpi_rank == 0) {
    hid_t                    file_id       = THDF_open_file(filename);
    THDF_3D_field_handler_t *field_handler = THDF_create_3D_field_handler_mpi(file_id,                            //
                                                                              normalized_output,                  //
                                                                              N_x, N_y, N_z,                      //
                                                                              N_incl, N_azimuth, N_frequencies);  //

    THDF_frequencies_grid_t frequencies_grid;
    frequencies_grid.N_frequencies = N_frequencies;
    frequencies_grid.frequencies   = frequencies;
    THDF_write_frequencies_grid_to_hdf5(file_id, &frequencies_grid);

    THDF_angular_grid_t angular_grid;
    angular_grid.N_azimuthal_angles   = N_azimuth;
    angular_grid.N_inclination_angles = N_incl;
    angular_grid.N_directions         = N_incl * N_azimuth;

    angular_grid.inclination_angles   = theta;
    angular_grid.inclinations_indices = inclinations_indices;
    angular_grid.azimuthal_angles     = chi;
    angular_grid.azimuthal_indices    = azimuthal_indices;
    THDF_write_angular_grid_to_hdf5(file_id, &angular_grid);

    THDF_close_3D_field_handler_mpi(field_handler);
    H5Fclose(file_id);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  int N_local_x, N_local_y, N_local_z;
  int local_start_x, local_start_y, local_start_z;
  int local_end_x, local_end_y, local_end_z;

  decompose_domain_3d(N_x, N_y, N_z, MPI_COMM_WORLD, &N_local_x, &N_local_y, &N_local_z, &local_start_x, &local_start_y,
                      &local_start_z, NULL, NULL, NULL);
  local_end_x = local_start_x + N_local_x;
  local_end_y = local_start_y + N_local_y;
  local_end_z = local_start_z + N_local_z;

  // if (mpi_rank >= 0) {
  printf("Rank %d: Local domain in X [%d:%d], Y [%d:%d], Z [%d:%d]\n", mpi_rank, local_start_x, local_end_x,
         local_start_y, local_end_y, local_start_z, local_end_z);
  // }

  int          stride_x, stride_y, stride_z, stride_incl, stride_azimuth, stride_frequencies, stride_stokes;
  const size_t input_size  = (size_t)N_local_x * N_local_y * N_local_z * N_incl * N_azimuth * N_frequencies * N_stokes;
  const size_t output_size = (size_t)N_local_x * N_local_y * N_local_z * N_incl * N_azimuth * N_frequencies;

  double         *stokes_IQUI = (double *)malloc(input_size * sizeof(double));
  THDF_float32_t *stokes_I    = (THDF_float32_t *)malloc(output_size * sizeof(THDF_float32_t));
  THDF_float32_t *stokes_QI   = (THDF_float32_t *)malloc(output_size * sizeof(THDF_float32_t));
  THDF_float32_t *stokes_UI   = (THDF_float32_t *)malloc(output_size * sizeof(THDF_float32_t));
  THDF_float32_t *stokes_VI   = (THDF_float32_t *)malloc(output_size * sizeof(THDF_float32_t));

  THDF_float_t *norm_multiplier_I  = NULL;
  THDF_float_t *norm_multiplier_QI = NULL;
  THDF_float_t *norm_multiplier_UI = NULL;
  THDF_float_t *norm_multiplier_VI = NULL;

  if (normalized_output) {
    const size_t norm_size = (size_t)N_local_x * N_local_y * N_local_z * N_incl * N_azimuth;
    norm_multiplier_I      = (THDF_float_t *)malloc(norm_size * sizeof(THDF_float_t));
    norm_multiplier_QI     = (THDF_float_t *)malloc(norm_size * sizeof(THDF_float_t));
    norm_multiplier_UI     = (THDF_float_t *)malloc(norm_size * sizeof(THDF_float_t));
    norm_multiplier_VI     = (THDF_float_t *)malloc(norm_size * sizeof(THDF_float_t));
  }

  make_input_example_dataset(stokes_IQUI, N_local_x, N_local_y, N_local_z, N_incl, N_azimuth, N_frequencies, N_stokes,
                             &stride_x, &stride_y, &stride_z, &stride_incl, &stride_azimuth, &stride_frequencies,
                             &stride_stokes);

  write_3d_field_block_mpi(filename,            //
                           MPI_COMM_WORLD,      //
                           normalized_output,   //
                           N_x,                 //
                           N_y,                 //
                           N_z,                 //
                           N_incl,              //
                           N_azimuth,           //
                           N_frequencies,       //
                           N_stokes,            //
                           N_local_x,           //
                           N_local_y,           //
                           N_local_z,           //
                           local_start_x,       //
                           local_start_y,       //
                           local_start_z,       //
                           stokes_IQUI,         //
                           stokes_I,            //
                           stokes_QI,           //
                           stokes_UI,           //
                           stokes_VI,           //
                           norm_multiplier_I,   //
                           norm_multiplier_QI,  //
                           norm_multiplier_UI,  //
                           norm_multiplier_VI,  //
                           stride_x,            //
                           stride_y,            //
                           stride_z,            //
                           stride_incl,         //
                           stride_azimuth,      //
                           stride_frequencies,  //
                           stride_stokes);      //

  free(stokes_IQUI);
  free(stokes_I);
  free(stokes_QI);
  free(stokes_UI);
  free(stokes_VI);
  free(norm_multiplier_I);
  free(norm_multiplier_QI);
  free(norm_multiplier_UI);
  free(norm_multiplier_VI);

  free(frequencies);
  free(theta);
  free(chi);
  free(inclinations_indices);
  free(azimuthal_indices);

  MPI_Finalize();
  return 0;
}

/////////////////////////////////////////////////////////////////////////
// main_3d_example_v2
/////////////////////////////////////////////////////////////////////////
int
main_3d_example_v2(int argc, char **argv) {

  (void)argc;
  (void)argv;

  MPI_Init(&argc, &argv);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  const int  N_x               = read_env_int_or_default("HDF_N_X", 8);
  const int  N_y               = read_env_int_or_default("HDF_N_Y", 8);
  const int  N_z               = read_env_int_or_default("HDF_N_Z", 17);
  const int  N_frequencies     = read_env_int_or_default("HDF_N_FREQUENCIES", 96);
  const int  N_incl            = read_env_int_or_default("HDF_N_INCL", 8);
  const int  N_azimuth         = read_env_int_or_default("HDF_N_AZIMUTH", 16);
  const int  N_stokes          = read_env_int_or_default("HDF_N_STOKES", 4);
  const bool normalized_output = read_env_int_or_default("HDF_NORMALIZE_OUTPUT", 0) != 0;
  int        step_z            = read_env_int_or_default("HDF_STEP_Z", 1);

  double *frequencies = (double *)malloc(N_frequencies * sizeof(double));

  for (int i = 0; i < N_frequencies; i++) {
    frequencies[i] = 1.0e14 + i * 1.0e12;
  }

  double *theta                = (double *)malloc(N_incl * sizeof(double));
  int    *inclinations_indices = (int *)malloc(N_incl * N_azimuth * sizeof(int));

  double *chi               = (double *)malloc(N_azimuth * sizeof(double));
  int    *azimuthal_indices = (int *)malloc(N_incl * N_azimuth * sizeof(int));

  double *heights = (double *)malloc(N_z * sizeof(double));
  for (int k = 0; k < N_z; k++) {
    heights[k] = k * 100.0;  // example heights in km
  }

  for (int i = 0; i < N_incl; i++) {
    theta[i]                = (3.141592653589793 / (N_incl - 1)) * (i);  // example values from 0 to 90 degrees in radians
    inclinations_indices[i] = i;
  }

  for (int j = 0; j < N_azimuth; j++) {
    chi[j]               = j * (3.141592653589793 * 2.0 / N_azimuth);  // example values from 0 to 360 degrees in radians
    azimuthal_indices[j] = j;
  }

  for (int i = 0; i < N_incl; i++) {
    for (int j = 0; j < N_azimuth; j++) {
      inclinations_indices[i * N_azimuth + j] = i;
      azimuthal_indices[i * N_azimuth + j]    = j;
    }
  }

  char filename[PATH_MAX];
  build_output_filepath(filename, sizeof(filename), "output_3D_field_mpi.h5");

  if (mpi_rank == 0) {
    hid_t                    file_id       = THDF_open_file(filename);
    THDF_3D_field_handler_t *field_handler = THDF_create_3D_field_handler_mpi(file_id,                            //
                                                                              normalized_output,                  //
                                                                              N_x, N_y, N_z,                      //
                                                                              N_incl, N_azimuth, N_frequencies);  //

    THDF_frequencies_grid_t frequencies_grid;
    frequencies_grid.N_frequencies = N_frequencies;
    frequencies_grid.frequencies   = frequencies;
    THDF_write_frequencies_grid_to_hdf5(file_id, &frequencies_grid);

    THDF_angular_grid_t angular_grid;
    angular_grid.N_azimuthal_angles   = N_azimuth;
    angular_grid.N_inclination_angles = N_incl;
    angular_grid.N_directions         = N_incl * N_azimuth;

    angular_grid.inclination_angles   = theta;
    angular_grid.inclinations_indices = inclinations_indices;
    angular_grid.azimuthal_angles     = chi;
    angular_grid.azimuthal_indices    = azimuthal_indices;
    THDF_write_angular_grid_to_hdf5(file_id, &angular_grid);

    THDF_geometry_3D_t geometry_3D;
    geometry_3D.N_x     = N_x;
    geometry_3D.N_y     = N_y;
    geometry_3D.N_z     = N_z;
    geometry_3D.heights = heights;
    geometry_3D.delta   = 100.0;
    THDF_write_geometry_3D_to_hdf5(file_id, &geometry_3D);

    THDF_close_3D_field_handler_mpi(field_handler);
    H5Fclose(file_id);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  int N_local_x, N_local_y, N_local_z;
  int local_start_x, local_start_y, local_start_z;
  int local_end_x, local_end_y, local_end_z;

  decompose_domain_3d(N_x, N_y, N_z, MPI_COMM_WORLD, &N_local_x, &N_local_y, &N_local_z, &local_start_x, &local_start_y,
                      &local_start_z, NULL, NULL, NULL);

  local_end_x = local_start_x + N_local_x;
  local_end_y = local_start_y + N_local_y;
  local_end_z = local_start_z + N_local_z;

  if (mpi_rank == 0) {
    printf("Rank %d: Local domain in X [%d:%d], Y [%d:%d], Z [%d:%d]\n", mpi_rank, local_start_x, local_end_x,
           local_start_y, local_end_y, local_start_z, local_end_z);
  }

  int          stride_x, stride_y, stride_z, stride_incl, stride_azimuth, stride_frequencies, stride_stokes;
  const size_t input_size  = (size_t)N_local_x * N_local_y * step_z * N_incl * N_azimuth * N_frequencies * N_stokes;
  const size_t output_size = (size_t)N_local_x * N_local_y * step_z * N_incl * N_azimuth * N_frequencies;

  double         *stokes_IQUI = (double *)malloc(input_size * sizeof(double));
  THDF_float32_t *stokes_I    = (THDF_float32_t *)malloc(output_size * sizeof(THDF_float32_t));
  THDF_float32_t *stokes_QI   = (THDF_float32_t *)malloc(output_size * sizeof(THDF_float32_t));
  THDF_float32_t *stokes_UI   = (THDF_float32_t *)malloc(output_size * sizeof(THDF_float32_t));
  THDF_float32_t *stokes_VI   = (THDF_float32_t *)malloc(output_size * sizeof(THDF_float32_t));

  THDF_float_t *norm_multiplier_I  = NULL;
  THDF_float_t *norm_multiplier_QI = NULL;
  THDF_float_t *norm_multiplier_UI = NULL;
  THDF_float_t *norm_multiplier_VI = NULL;

  if (normalized_output) {
    const size_t norm_size = (size_t)N_local_x * N_local_y * step_z * N_incl * N_azimuth;
    norm_multiplier_I      = (THDF_float_t *)malloc(norm_size * sizeof(THDF_float_t));
    norm_multiplier_QI     = (THDF_float_t *)malloc(norm_size * sizeof(THDF_float_t));
    norm_multiplier_UI     = (THDF_float_t *)malloc(norm_size * sizeof(THDF_float_t));
    norm_multiplier_VI     = (THDF_float_t *)malloc(norm_size * sizeof(THDF_float_t));
  }

  int max_N_local_z;
  int min_N_local_z;
  MPI_Allreduce(&N_local_z, &max_N_local_z, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&N_local_z, &min_N_local_z, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  if (max_N_local_z < step_z) {
    if (mpi_rank == 0) {
      printf("Warning: max local Z size (%d) is smaller than step_z (%d). Adjusting step_z to %d.\n", max_N_local_z,
             step_z, max_N_local_z);
    }
    step_z = 1;
  }

  if (mpi_rank == 0) {
    printf("Max local Z size across all ranks: %d\n", max_N_local_z);
    printf("Min local Z size across all ranks: %d\n", min_N_local_z);
    printf("Writing in Z blocks of size %d (last block may be smaller)\n", step_z);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int local_zi = 0; local_zi < max_N_local_z; local_zi += 1) {

    const int current_N_local_z = 1;
    const int is_writer         = (local_zi < N_local_z) ? 1 : 0;
    // Only ranks that still have Z slices to write in this block participate as writers

    const int id_z = local_start_z + local_zi;

    MPI_Comm write_communicator;
    MPI_Comm_split(MPI_COMM_WORLD, is_writer ? 1 : MPI_UNDEFINED, mpi_rank, &write_communicator);

    if (is_writer) {

      if (mpi_rank < 10) {
        printf("Rank %d writing Z block [%d:%d] (local Z indices [%d:%d])\n", mpi_rank, local_start_z + local_zi,
               local_start_z + local_zi + current_N_local_z, local_zi, local_zi + current_N_local_z);
      }

      make_input_example_dataset_zz(stokes_IQUI, N_local_x, N_local_y, current_N_local_z, N_incl, N_azimuth,
                                    N_frequencies, N_stokes, &stride_x, &stride_y, &stride_z, &stride_incl,
                                    &stride_azimuth, &stride_frequencies, &stride_stokes, id_z);

      write_3d_field_block_mpi(filename,            //
                               write_communicator,  //
                               normalized_output,   //
                               N_x,                 //
                               N_y,                 //
                               N_z,                 //
                               N_incl,              //
                               N_azimuth,           //
                               N_frequencies,       //
                               N_stokes,            //
                               N_local_x,           //
                               N_local_y,           //
                               1,                   //
                               local_start_x,       //
                               local_start_y,       //
                               id_z,                //
                               stokes_IQUI,         //
                               stokes_I,            //
                               stokes_QI,           //
                               stokes_UI,           //
                               stokes_VI,           //
                               norm_multiplier_I,   //
                               norm_multiplier_QI,  //
                               norm_multiplier_UI,  //
                               norm_multiplier_VI,  //
                               stride_x,            //
                               stride_y,            //
                               stride_z,            //
                               stride_incl,         //
                               stride_azimuth,      //
                               stride_frequencies,  //
                               stride_stokes);      //

      MPI_Barrier(write_communicator);  // Ensure all writers have finished before we free the communicator
      MPI_Comm_free(&write_communicator);
    }

    // Keep all ranks in lockstep between Z blocks; non-writers have MPI_COMM_NULL.
    MPI_Barrier(MPI_COMM_WORLD);
  }

  free(stokes_IQUI);
  free(stokes_I);
  free(stokes_QI);
  free(stokes_UI);
  free(stokes_VI);
  free(norm_multiplier_I);
  free(norm_multiplier_QI);
  free(norm_multiplier_UI);
  free(norm_multiplier_VI);
  free(heights);
  free(frequencies);
  free(theta);
  free(chi);
  free(inclinations_indices);
  free(azimuthal_indices);

  MPI_Finalize();
  return 0;
}  // END: main_3d_example_v2
