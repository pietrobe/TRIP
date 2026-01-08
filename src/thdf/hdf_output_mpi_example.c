#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "thdf.h"

/////////////////////////////////////////////////////////////
/// make_example_output_field_surface
/////////////////////////////////////////////////////////////
THDF_field_t
make_example_output_field_surface(int N_frequencies, const int N_incl, const int N_azimuth) {
  // Create output field for a single (x, y) position with all inclinations and azimuths
  // The Stokes arrays must be contiguous: size = N_incl * N_azimuth * N_frequencies
  const int angular_size = N_incl * N_azimuth;

  THDF_field_t output_field;
  output_field.index_i       = 0;
  output_field.index_j       = 0;
  output_field.index_incl    = 0;
  output_field.index_azimuth = 0;

  // Allocate contiguous arrays for all angular directions and frequencies
  output_field.stokes_I  = (THDF_float_t *)malloc(angular_size * N_frequencies * sizeof(THDF_float_t));
  output_field.stokes_QI = (THDF_float_t *)malloc(angular_size * N_frequencies * sizeof(THDF_float_t));
  output_field.stokes_UI = (THDF_float_t *)malloc(angular_size * N_frequencies * sizeof(THDF_float_t));
  output_field.stokes_VI = (THDF_float_t *)malloc(angular_size * N_frequencies * sizeof(THDF_float_t));

  // Fill with example data: layout is [incl][azimuth][frequencies]
  for (int incl = 0; incl < N_incl; incl++) {
    for (int azim = 0; azim < N_azimuth; azim++) {
      int angular_idx = incl * N_azimuth + azim;
      for (int f = 0; f < N_frequencies; f++) {
        int idx                     = angular_idx * N_frequencies + f;
        output_field.stokes_I[idx]  = 1.0 * (f + 1) * (incl + 1) * (azim + 1) * 0.01;
        output_field.stokes_QI[idx] = 2.0 * (f + 1) * (incl + 1) * (azim + 1) * 0.01;
        output_field.stokes_UI[idx] = 3.0 * (f + 1) * (incl + 1) * (azim + 1) * 0.01;
        output_field.stokes_VI[idx] = 4.0 * (f + 1) * (incl + 1) * (azim + 1) * 0.01;
      }
    }
  }

  return output_field;
}

/////////////////////////////////////////////////////////////
/// make_example_output_field_surface_sd
/////////////////////////////////////////////////////////////
THDF_field_t
make_example_output_field_surface_sd(int N_frequencies, const int N_incl, const int N_azimuth, const int N_x,
                                     const int N_y) {
  // Create output field for a single (x, y) position with all inclinations and azimuths
  // The Stokes arrays must be contiguous: size = N_incl * N_azimuth * N_frequencies
  const int angular_size = N_incl * N_azimuth;

  THDF_field_t output_field;
  output_field.index_i       = 0;
  output_field.index_j       = 0;
  output_field.index_incl    = 0;
  output_field.index_azimuth = 0;

  // Allocate contiguous arrays for all angular directions and frequencies
  output_field.stokes_I  = (THDF_float_t *)malloc(angular_size * N_frequencies * N_x * N_y * sizeof(THDF_float_t));
  output_field.stokes_QI = (THDF_float_t *)malloc(angular_size * N_frequencies * N_x * N_y * sizeof(THDF_float_t));
  output_field.stokes_UI = (THDF_float_t *)malloc(angular_size * N_frequencies * N_x * N_y * sizeof(THDF_float_t));
  output_field.stokes_VI = (THDF_float_t *)malloc(angular_size * N_frequencies * N_x * N_y * sizeof(THDF_float_t));

  // Fill with example data: layout is [incl][azimuth][frequencies]
  for (int x = 0; x < N_x; x++) {
    for (int y = 0; y < N_y; y++) {
      for (int incl = 0; incl < N_incl; incl++) {
        for (int azim = 0; azim < N_azimuth; azim++) {
          int angular_idx = incl * N_azimuth + azim;
          for (int f = 0; f < N_frequencies; f++) {
            int idx                     = ((x * N_y + y) * angular_size + angular_idx) * N_frequencies + f;
            output_field.stokes_I[idx]  = sin(1.0 * (f + 1) * (incl + 1) * (azim + 1) * 0.5 * (x + 1) * (y + 1));
            output_field.stokes_QI[idx] = sin(2.0 * (f + 1) * (incl + 1) * (azim + 1) * 0.01 * (x + 1) * (y + 1));
            output_field.stokes_UI[idx] = cos(3.0 * (f + 1) * (incl + 1) * (azim + 1) * 0.01 * (x + 1) * (y + 1));
            output_field.stokes_VI[idx] = cos(1.0 * (f + 1) * (incl + 1) * (azim + 1) * -0.05 * (x + 1) * (y + 1));
          }
        }
      }
    }
  }

  return output_field;
}

/////////////////////////////////////////////////////////////
/// destroy_example_output_field_surface
/////////////////////////////////////////////////////////////
void
destroy_example_output_field_surface(THDF_field_t *output_field) {
  if (output_field == NULL) {
    return;
  }

  free(output_field->stokes_I);
  free(output_field->stokes_QI);
  free(output_field->stokes_UI);
  free(output_field->stokes_VI);
}

//============================================================================
// main
//============================================================================
/////////////////////////////////////////////////////////////
/// main_default
/////////////////////////////////////////////////////////////
int
main_default(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int         N_x;
  int         N_y;
  int         N_z;
  int         N_theta;
  int         N_chi;
  int         N_frequencies;
  const char *env;

  /* defaults */
  N_x = 3;
  env = getenv("HDF5_N_X");
  if (env)
    N_x = atoi(env);

  N_y = 3;
  env = getenv("HDF5_N_Y");
  if (env)
    N_y = atoi(env);

  N_z = 2;
  env = getenv("HDF5_N_Z");
  if (env)
    N_z = atoi(env);

  N_theta = 6;
  env     = getenv("HDF5_N_THETA");
  if (env)
    N_theta = atoi(env);

  N_chi = 8;
  env   = getenv("HDF5_N_CHI");
  if (env)
    N_chi = atoi(env);

  N_frequencies = 130;
  env           = getenv("HDF5_N_FREQUENCIES");
  if (env)
    N_frequencies = atoi(env);
  double *frequencies = NULL;

  const char filename[100] = "output_field_mpi.h5";

  int mpi_index_x = mpi_rank % N_x;
  int mpi_index_y = (mpi_rank / N_x) % N_y;
  int mpi_index_z = (mpi_rank / (N_x * N_y)) % N_z;
  //   int mpi_surface_size = N_x * N_y;

  // Each rank creates its own local output field with contiguous Stokes arrays
  THDF_field_t output_field = make_example_output_field_surface(N_frequencies, N_theta, N_chi);

  for (int rank = 0; rank < mpi_size; rank++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == rank) {
      printf("Rank %d/%d: mpi_index_x=%d, mpi_index_y=%d, mpi_index_z=%d\n", mpi_rank, mpi_size, mpi_index_x, mpi_index_y,
             mpi_index_z);
      fflush(stdout);
    }
  }

  frequencies = (double *)malloc(N_frequencies * sizeof(double));
  for (int i = 0; i < N_frequencies; i++) {
    frequencies[i] = 5.0e14 + i * 1.0e11;  // Example: frequencies in Hz
  }

  double *theta                = NULL;
  double *chi                  = NULL;
  int    *inclinations_indices = NULL;
  int    *azimuthal_indices    = NULL;

  theta                = (double *)malloc(N_theta * sizeof(double));
  chi                  = (double *)malloc(N_chi * sizeof(double));
  inclinations_indices = (int *)malloc(N_theta * sizeof(int));
  azimuthal_indices    = (int *)malloc(N_chi * sizeof(int));

  for (int i = 0; i < N_theta; i++) {
    theta[i]                = i * 10.0;  // just example values
    inclinations_indices[i] = i;
  }

  for (int j = 0; j < N_chi; j++) {
    chi[j]               = j * 15.0;  // just example values
    azimuthal_indices[j] = j;
  }

  if (mpi_rank == 0) {
    printf("Creating MPI-enabled output field dataset... from rank 0\n");

    hid_t file_id, plist_id;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id  = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    // write frequencies grid
    THDF_frequencies_grid_t freq_grid;
    freq_grid.N_frequencies = N_frequencies;
    freq_grid.frequencies   = frequencies;
    if (THDF_write_frequencies_grid_to_hdf5(file_id, &freq_grid) != 0) {
      fprintf(stderr, "Error writing frequencies grid to HDF5 file\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    THDF_angular_grid_t angular_grid;
    angular_grid.N_directions         = N_theta * N_chi;
    angular_grid.N_inclination_angles = N_theta;
    angular_grid.N_azimuthal_angles   = N_chi;
    angular_grid.inclination_angles   = theta;
    angular_grid.azimuthal_angles     = chi;
    angular_grid.inclinations_indices = inclinations_indices;
    angular_grid.azimuthal_indices    = azimuthal_indices;

    if (THDF_write_angular_grid_to_hdf5(file_id, &angular_grid) != 0) {
      fprintf(stderr, "Error writing emergent angular grid to HDF5 file\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    H5Fclose(file_id);
  }

  // All ranks must synchronize before parallel file open
  MPI_Barrier(MPI_COMM_WORLD);

  // Create a sub-communicator for ranks that will write (mpi_index_z == 0)
  int      is_writer = (mpi_index_z == 0) ? 1 : 0;
  MPI_Comm writer_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_writer, mpi_rank, &writer_comm);

  // Only the writing ranks participate in parallel file operations
  if (mpi_index_z == 0) {
    printf("Rank %d: Writing to MPI-enabled output field dataset... to index (%d, %d)\n", mpi_rank, mpi_index_x,
           mpi_index_y);
    fflush(stdout);

    hid_t file_id, plist_id;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, writer_comm, MPI_INFO_NULL);
    file_id = H5Fopen(filename, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    THDF_field_handler_t *output_dset_handler =
        THDF_create_field_handler_mpi(file_id, N_x, N_y, N_theta, N_chi, N_frequencies);

    if (output_dset_handler == NULL) {
      fprintf(stderr, "Rank %d: Error creating MPI output field dataset\n", mpi_rank);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // Write local output field surface data
    THDF_write_field_dataset_to_hdf5(output_dset_handler,                         //
                                     &output_field,                               //
                                     mpi_index_x,                                 //
                                     mpi_index_y,                                 //
                                     0, 0, 1, 1, N_theta, N_chi, N_frequencies);  //

    THDF_close_field_handler_mpi(output_dset_handler);
    H5Fclose(file_id);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_free(&writer_comm);

  MPI_Finalize();

  free(frequencies);
  free(theta);
  free(chi);
  free(inclinations_indices);
  free(azimuthal_indices);
  destroy_example_output_field_surface(&output_field);

  return 0;
}

THDF_JKQ_field_t *
make_example_JKQ_field(const int N_x_local, const int N_y_local, const int N_z_local, const int N_JKQ_real,
                       const int N_JKQ_imag, const int N_frequencies) {

  THDF_JKQ_field_t *JKQ_field = (THDF_JKQ_field_t *)malloc(sizeof(THDF_JKQ_field_t));

  if (!JKQ_field) {
    fprintf(stderr, "Error allocating memory for JKQ field\n");
    return NULL;
  }

  const int total_size_real      = N_x_local * N_y_local * N_z_local * N_JKQ_real * N_frequencies;
  const int total_size_imag      = N_x_local * N_y_local * N_z_local * N_JKQ_imag * N_frequencies;
  const int total_norm_size_real = N_x_local * N_y_local * N_z_local * N_JKQ_real;
  const int total_norm_size_imag = N_x_local * N_y_local * N_z_local * N_JKQ_imag;

  JKQ_field->JKQ_re                 = (THDF_JKQ_float_t *)malloc(total_size_real * sizeof(THDF_JKQ_float_t));
  JKQ_field->JKQ_im                 = (THDF_JKQ_float_t *)malloc(total_size_imag * sizeof(THDF_JKQ_float_t));
  JKQ_field->norm_multiplier_JKQ_re = (THDF_JKQ_n_float_t *)malloc(total_norm_size_real * sizeof(THDF_JKQ_n_float_t));
  JKQ_field->norm_multiplier_JKQ_im = (THDF_JKQ_n_float_t *)malloc(total_norm_size_imag * sizeof(THDF_JKQ_n_float_t));

  if (!JKQ_field->JKQ_re || !JKQ_field->JKQ_im || !JKQ_field->norm_multiplier_JKQ_re ||
      !JKQ_field->norm_multiplier_JKQ_im) {
    fprintf(stderr, "Error allocating memory for JKQ field arrays\n");
    free(JKQ_field->JKQ_re);
    free(JKQ_field->JKQ_im);
    free(JKQ_field->norm_multiplier_JKQ_re);
    free(JKQ_field->norm_multiplier_JKQ_im);
    free(JKQ_field);
    return NULL;
  }

  // Fill with example data
  for (int i = 0; i < total_size_real; i++) {
    JKQ_field->JKQ_re[i] = (THDF_JKQ_float_t)(0.1 * (i + 1));
  }
  for (int i = 0; i < total_size_imag; i++) {
    JKQ_field->JKQ_im[i] = (THDF_JKQ_float_t)(-0.1 * (i + 1));
  }
  for (int i = 0; i < total_norm_size_real; i++) {
    JKQ_field->norm_multiplier_JKQ_re[i] = (THDF_JKQ_n_float_t)(1.0 + 0.01 * i);
  }
  for (int i = 0; i < total_norm_size_imag; i++) {
    JKQ_field->norm_multiplier_JKQ_im[i] = (THDF_JKQ_n_float_t)(1.0 - 0.01 * i);
  }

  return JKQ_field;
}

/////////////////////////////////////////////////////////////
/// destroy_example_JKQ_field
/////////////////////////////////////////////////////////////
void
destroy_example_JKQ_field(THDF_JKQ_field_t *JKQ_field) {
  if (JKQ_field == NULL) {
    return;
  }

  free(JKQ_field->JKQ_re);
  free(JKQ_field->JKQ_im);
  free(JKQ_field->norm_multiplier_JKQ_re);
  free(JKQ_field->norm_multiplier_JKQ_im);
  free(JKQ_field);
}

/////////////////////////////////////////////////////////////
/// main_domain_dec
/////////////////////////////////////////////////////////////
int
main_JKQ_domain_dec(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int   N_x;
  int   N_y;
  int   N_z;
  int   N_JKQ_real;
  int   N_JKQ_imag;
  int   N_frequencies;
  char *env;

  /* defaults */
  N_x = 15;
  env = getenv("HDF5_N_X");
  if (env)
    N_x = atoi(env);

  N_y = 15;
  env = getenv("HDF5_N_Y");
  if (env)
    N_y = atoi(env);

  N_z = 2;
  env = getenv("HDF5_N_Z");
  if (env)
    N_z = atoi(env);

  N_JKQ_real = 6;
  env        = getenv("HDF5_N_JKQ_REAL");
  if (env)
    N_JKQ_real = atoi(env);

  N_JKQ_imag = 3;
  env        = getenv("HDF5_N_JKQ_IMAG");
  if (env)
    N_JKQ_imag = atoi(env);

  N_frequencies = 96;
  env           = getenv("HDF5_N_FREQUENCIES");
  if (env)
    N_frequencies = atoi(env);

  double *frequencies                = NULL;
  int    *KQ_values                  = NULL;
  int    *KQ_values_commpressed_real = NULL;
  int    *KQ_values_commpressed_imag = NULL;

  // Domain decomposition setup
  int       ranks_per_surface = mpi_size;
  const int surface_rank      = mpi_rank % ranks_per_surface;  // position within surface
  // const int mpi_index_z       = mpi_rank / ranks_per_surface;  // which z-layer
  // const int is_surface_writer = (mpi_index_z == 0) ? 1 : 0;

  const int subdiv_x = (int)ceil(sqrt((double)ranks_per_surface));
  const int subdiv_y = (int)ceil((double)ranks_per_surface / subdiv_x);

  if (mpi_rank == 0)
    printf("Ranks per surface: %d, subdiv_x: %d, subdiv_y: %d, total domains on surface: %d\n", ranks_per_surface,
           subdiv_x, subdiv_y, subdiv_x * subdiv_y);

  if (subdiv_x * subdiv_y != ranks_per_surface) {
    ranks_per_surface = subdiv_x * subdiv_y;  // adjust to actual number of ranks used on surface
    if (mpi_rank == 0)
      printf("Adjusted ranks per surface to %d\n", ranks_per_surface);
  }

  int index_x = surface_rank % subdiv_x;
  int index_y = surface_rank / subdiv_x;

  const int base_N_x = N_x / subdiv_x;
  const int rem_x    = N_x % subdiv_x;
  const int base_N_y = N_y / subdiv_y;
  const int rem_y    = N_y % subdiv_y;

  int local_N_x     = base_N_x + (index_x < rem_x ? 1 : 0);
  int local_N_y     = base_N_y + (index_y < rem_y ? 1 : 0);
  int local_start_x = index_x * base_N_x + (index_x < rem_x ? index_x : rem_x);
  int local_start_y = index_y * base_N_y + (index_y < rem_y ? index_y : rem_y);

  const char filename[100] = "output_JKQ_mpi.h5";

  frequencies = (double *)malloc(N_frequencies * sizeof(double));
  for (int i = 0; i < N_frequencies; i++) {
    frequencies[i] = 5.0e14 + i * 1.0e11;  // Example: frequencies in Hz
  }

  KQ_values = (int *)malloc(9 * 2 * sizeof(int));  //
  // Fill KQ_values with example data
  for (int i = 0; i < 9 * 2; i++) KQ_values[i] = i + 1;

  KQ_values_commpressed_real = (int *)malloc(N_JKQ_real * 2 * sizeof(int));  //
  KQ_values_commpressed_imag = (int *)malloc(N_JKQ_imag * 2 * sizeof(int));  //
  // Fill compressed KQ values with example data
  for (int i = 0; i < N_JKQ_real * 2; i++) KQ_values_commpressed_real[i] = i + 1;
  for (int i = 0; i < N_JKQ_imag * 2; i++) KQ_values_commpressed_imag[i] = i + 11;

  if (mpi_rank == 0) {
    printf("Creating MPI-enabled JKQ dataset... from rank 0\n");

    hid_t file_id, plist_id;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id  = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    // write frequencies grid
    THDF_frequencies_grid_t freq_grid;
    freq_grid.N_frequencies = N_frequencies;
    freq_grid.frequencies   = frequencies;
    if (THDF_write_frequencies_grid_to_hdf5(file_id, &freq_grid) != 0) {
      fprintf(stderr, "Error writing frequencies grid to HDF5 file\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // write JKQ KQ matrix
    THDF_KQ_matrix_t KQ_matrix;
    KQ_matrix.KQ_size                 = 9;
    KQ_matrix.KQ_matrix               = KQ_values;
    KQ_matrix.KQ_compressed_size_real = N_JKQ_real;
    KQ_matrix.KQ_compressed_size_imag = N_JKQ_imag;
    KQ_matrix.KQ_compressed_real      = KQ_values_commpressed_real;
    KQ_matrix.KQ_compressed_imag      = KQ_values_commpressed_imag;

    if (THDF_write_KQ_matrix_to_hdf5(file_id, &KQ_matrix) != 0) {
      fprintf(stderr, "Error writing KQ matrix to HDF5 file\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    H5Fclose(file_id);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  //// free datasets and communicators
  free(frequencies);
  free(KQ_values);
  free(KQ_values_commpressed_real);
  free(KQ_values_commpressed_imag);

  {
    printf("Rank %2d: Writing to MPI-enabled output JKQ dataset... from index (%2d, %2d, %2d) to end (%2d, %2d, %2d)\n",
           mpi_rank, local_start_x, local_start_y, 0, local_start_x + local_N_x - 1, local_start_y + local_N_y - 1,
           N_z - 1);
    fflush(stdout);

    hid_t file_id, plist_id;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    file_id = H5Fopen(filename, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    THDF_JKQ_handler_t *JKQ_dset_handler =
        THDF_create_JKQ_field_handler_mpi(file_id, N_x, N_y, N_z, N_JKQ_real, N_JKQ_imag, N_frequencies);

    if (JKQ_dset_handler == NULL) {
      fprintf(stderr, "Rank %d: Error creating MPI JKQ dataset\n", mpi_rank);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    THDF_JKQ_field_t *JKQ_field =
        make_example_JKQ_field(local_N_x, local_N_y, N_z, N_JKQ_real, N_JKQ_imag, N_frequencies);

    // Write local JKQ field data
    THDF_write_JKQ_field_to_hdf5(JKQ_dset_handler,  //
                                 JKQ_field,         //
                                 local_start_x,     //
                                 local_start_y,     //
                                 0,                 //
                                 local_N_x,         //
                                 local_N_y,         //
                                 N_z);              //

    THDF_close_JKQ_field_handler_mpi(JKQ_dset_handler);
    H5Fclose(file_id);

    destroy_example_JKQ_field(JKQ_field);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}

/////////////////////////////////////////////////////////////
/// main_domain_dec
/////////////////////////////////////////////////////////////
int
main_domain_dec(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int   N_x;
  int   N_y;
  int   N_z;
  int   N_theta;
  int   N_chi;
  int   N_frequencies;
  char *env;

  /* defaults */
  N_x = 15;
  env = getenv("HDF5_N_X");
  if (env)
    N_x = atoi(env);

  N_y = 15;
  env = getenv("HDF5_N_Y");
  if (env)
    N_y = atoi(env);

  N_z = 2;
  env = getenv("HDF5_N_Z");
  if (env)
    N_z = atoi(env);

  N_theta = 6;
  env     = getenv("HDF5_N_THETA");
  if (env)
    N_theta = atoi(env);

  N_chi = 8;
  env   = getenv("HDF5_N_CHI");
  if (env)
    N_chi = atoi(env);

  N_frequencies = 130;
  env           = getenv("HDF5_N_FREQUENCIES");
  if (env)
    N_frequencies = atoi(env);
  double *frequencies = NULL;

  const char filename[100] = "output_field_mpi.h5";

  // const int mpi_surface_size = N_x * N_y;

  const int ranks_per_surface = mpi_size / N_z;
  const int surface_rank      = mpi_rank % ranks_per_surface;  // position within surface
  const int mpi_index_z       = mpi_rank / ranks_per_surface;  // which z-layer
  const int is_surface_writer = (mpi_index_z == 0) ? 1 : 0;

  const int subdiv_x = (int)ceil(sqrt((double)ranks_per_surface));
  const int subdiv_y = (int)ceil((double)ranks_per_surface / subdiv_x);

  int index_x = surface_rank % subdiv_x;
  int index_y = surface_rank / subdiv_x;

  int local_N_x     = N_x / subdiv_x;
  int local_N_y     = N_y / subdiv_y;
  int local_start_x = index_x * local_N_x;
  int local_start_y = index_y * local_N_y;

  if (index_x == subdiv_x - 1) {
    local_N_x = N_x - local_start_x;  // last rank takes the remainder
  }

  if (index_y == subdiv_y - 1) {
    local_N_y = N_y - local_start_y;  // last rank takes the remainder
  }

  frequencies = (double *)malloc(N_frequencies * sizeof(double));
  for (int i = 0; i < N_frequencies; i++) {
    frequencies[i] = 5.0e14 + i * 1.0e11;  // Example: frequencies in Hz
  }

  double *theta                = NULL;
  double *chi                  = NULL;
  int    *inclinations_indices = NULL;
  int    *azimuthal_indices    = NULL;

  theta                = (double *)malloc(N_theta * sizeof(double));
  chi                  = (double *)malloc(N_chi * sizeof(double));
  inclinations_indices = (int *)malloc(N_theta * sizeof(int));
  azimuthal_indices    = (int *)malloc(N_chi * sizeof(int));

  for (int i = 0; i < N_theta; i++) {
    theta[i]                = i * 10.0;  // just example values
    inclinations_indices[i] = i;
  }

  for (int j = 0; j < N_chi; j++) {
    chi[j]               = j * 15.0;  // just example values
    azimuthal_indices[j] = j;
  }

  if (mpi_rank == 0) {
    printf("Creating MPI-enabled output field dataset... from rank 0\n");

    hid_t file_id, plist_id;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id  = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    // write frequencies grid
    THDF_frequencies_grid_t freq_grid;
    freq_grid.N_frequencies = N_frequencies;
    freq_grid.frequencies   = frequencies;
    if (THDF_write_frequencies_grid_to_hdf5(file_id, &freq_grid) != 0) {
      fprintf(stderr, "Error writing frequencies grid to HDF5 file\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    THDF_angular_grid_t angular_grid;
    angular_grid.N_directions         = N_theta * N_chi;
    angular_grid.N_inclination_angles = N_theta;
    angular_grid.N_azimuthal_angles   = N_chi;
    angular_grid.inclination_angles   = theta;
    angular_grid.azimuthal_angles     = chi;
    angular_grid.inclinations_indices = inclinations_indices;
    angular_grid.azimuthal_indices    = azimuthal_indices;

    if (THDF_write_angular_grid_to_hdf5(file_id, &angular_grid) != 0) {
      fprintf(stderr, "Error writing emergent angular grid to HDF5 file\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    H5Fclose(file_id);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Create a sub-communicator for ranks that will write (mpi_index_z == 0)
  MPI_Comm writer_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_surface_writer, mpi_rank, &writer_comm);

  // Only the writing ranks participate in parallel file operations
  if (is_surface_writer) {
    printf("Rank %d: Writing to MPI-enabled output field dataset... from index (%d, %d) to end (%d, %d)\n", mpi_rank,
           local_start_x, local_start_y, local_start_x + local_N_x - 1, local_start_y + local_N_y - 1);
    fflush(stdout);

    hid_t file_id, plist_id;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, writer_comm, MPI_INFO_NULL);
    file_id = H5Fopen(filename, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    THDF_field_handler_t *output_dset_handler =
        THDF_create_field_handler_mpi(file_id, N_x, N_y, N_theta, N_chi, N_frequencies);
    if (output_dset_handler == NULL) {
      fprintf(stderr, "Rank %d: Error creating MPI output field dataset\n", mpi_rank);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // Each rank creates its own local output field with contiguous Stokes arrays
    THDF_field_t output_field = make_example_output_field_surface_sd(N_frequencies, N_theta, N_chi, local_N_x, local_N_y);
    // Write local output field surface data
    THDF_write_field_dataset_to_hdf5(output_dset_handler,                                         //
                                     &output_field,                                               //
                                     local_start_x,                                               //
                                     local_start_y,                                               //
                                     0, 0, local_N_x, local_N_y, N_theta, N_chi, N_frequencies);  //

    destroy_example_output_field_surface(&output_field);
    THDF_close_field_handler_mpi(output_dset_handler);
    H5Fclose(file_id);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_free(&writer_comm);
  MPI_Finalize();
  free(frequencies);
  free(theta);
  free(chi);
  free(inclinations_indices);
  free(azimuthal_indices);

  return 0;
}

/////////////////////////////////////////////////////////////
/// main_example_ar_domain_dec
/////////////////////////////////////////////////////////////
int
main_example_ar_domain_dec(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int   N_x;
  int   N_y;
  int   N_z;
  int   N_frequencies;
  int   N_DIRECTIONS;
  char *env;

  /* defaults */
  N_x = 15;
  env = getenv("HDF5_N_X");
  if (env)
    N_x = atoi(env);

  N_y = 15;
  env = getenv("HDF5_N_Y");
  if (env)
    N_y = atoi(env);

  N_z = 2;
  env = getenv("HDF5_N_Z");
  if (env)
    N_z = atoi(env);

  N_frequencies = 55;
  env           = getenv("HDF5_N_FREQUENCIES");
  if (env)
    N_frequencies = atoi(env);

  N_DIRECTIONS = 5;
  env          = getenv("HDF5_N_DIRECTIONS");
  if (env)
    N_DIRECTIONS = atoi(env);

  double *frequencies = NULL;
  double *muv         = NULL;
  double *chiv        = NULL;

  const char filename[100] = "output_field_ar_mpi.h5";

  // Domain decomposition setup
  int       ranks_per_surface = mpi_size / N_z;
  const int surface_rank      = mpi_rank % ranks_per_surface;  // position within surface
  const int mpi_index_z       = mpi_rank / ranks_per_surface;  // which z-layer
  const int is_surface_writer = (mpi_index_z == 0) ? 1 : 0;

  const int subdiv_x = (int)ceil(sqrt((double)ranks_per_surface));
  const int subdiv_y = (int)ceil((double)ranks_per_surface / subdiv_x);

  if (mpi_rank == 0)
    printf("Ranks per surface: %d, subdiv_x: %d, subdiv_y: %d, total domains on surface: %d\n", ranks_per_surface,
           subdiv_x, subdiv_y, subdiv_x * subdiv_y);

  if (subdiv_x * subdiv_y != ranks_per_surface) {
    ranks_per_surface = subdiv_x * subdiv_y;  // adjust to actual number of ranks used on surface
    if (mpi_rank == 0)
      printf("Adjusted ranks per surface to %d\n", ranks_per_surface);
  }

  int index_x = surface_rank % subdiv_x;
  int index_y = surface_rank / subdiv_x;

  const int base_N_x = N_x / subdiv_x;
  const int rem_x    = N_x % subdiv_x;
  const int base_N_y = N_y / subdiv_y;
  const int rem_y    = N_y % subdiv_y;

  int local_N_x     = base_N_x + (index_x < rem_x ? 1 : 0);
  int local_N_y     = base_N_y + (index_y < rem_y ? 1 : 0);
  int local_start_x = index_x * base_N_x + (index_x < rem_x ? index_x : rem_x);
  int local_start_y = index_y * base_N_y + (index_y < rem_y ? index_y : rem_y);

  // Allocate and initialize frequencies
  frequencies = (double *)malloc(N_frequencies * sizeof(double));
  for (int i = 0; i < N_frequencies; i++) {
    frequencies[i] = 1.0e14 + i * 1.0e12;
  }

  // Allocate and initialize arbitrary directions
  muv  = (double *)malloc(N_DIRECTIONS * sizeof(double));
  chiv = (double *)malloc(N_DIRECTIONS * sizeof(double));
  for (int i = 0; i < N_DIRECTIONS; i++) {
    muv[i]  = 0.1 + i * 0.1;
    chiv[i] = 0.0;
  }

  // Rank 0 creates the file and writes metadata
  if (mpi_rank == 0) {
    printf("Creating MPI-enabled ARD output field dataset... from rank 0\n");

    hid_t file_id, plist_id;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id  = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    // write frequencies grid
    THDF_ard_frequencies_t ard_freqs;
    ard_freqs.N_frequencies = N_frequencies;
    ard_freqs.frequencies   = frequencies;
    if (THDF_write_ard_frequencies(file_id, &ard_freqs) != 0) {
      fprintf(stderr, "Error writing ARD frequencies grid to HDF5 file\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // write arbitrary directions
    THDF_ard_directions_t ard_directions;
    ard_directions.N_directions = N_DIRECTIONS;
    ard_directions.mu           = muv;
    ard_directions.chi          = chiv;
    if (THDF_write_ard_directions(file_id, &ard_directions) != 0) {
      fprintf(stderr, "Error writing ARD directions to HDF5 file\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    H5Fclose(file_id);
  }

  // All ranks must synchronize before parallel file open
  MPI_Barrier(MPI_COMM_WORLD);

  // Create a sub-communicator for ranks that will write (mpi_index_z == 0)
  MPI_Comm writer_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_surface_writer, mpi_rank, &writer_comm);

  // Only the writing ranks participate in parallel file operations
  if (is_surface_writer) {
    printf("Rank %d: Writing to MPI-enabled ARD output field dataset... from index (%d, %d) to end (%d, %d)\n", mpi_rank,
           local_start_x, local_start_y, local_start_x + local_N_x - 1, local_start_y + local_N_y - 1);
    fflush(stdout);

    hid_t file_id, plist_id;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, writer_comm, MPI_INFO_NULL);
    file_id = H5Fopen(filename, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    // Loop over all directions
    for (int di = 0; di < N_DIRECTIONS; di++) {

      THDF_ard_field_t field_data;

      // Create ARD field handler for this direction
      THDF_ard_field_handler *dset_handler = THDF_create_ard_field_handler_mpi(file_id,       //
                                                                               muv[di],       //
                                                                               chiv[di],      //
                                                                               N_x,           //
                                                                               N_y,           //
                                                                               N_DIRECTIONS,  //
                                                                               N_frequencies);

      if (dset_handler == NULL) {
        fprintf(stderr, "Rank %d: Error creating ARD field handler for direction %d\n", mpi_rank, di);
        MPI_Abort(MPI_COMM_WORLD, -1);
      }

      // Loop over local subdomain
      for (int ix = 0; ix < local_N_x; ix++) {
        for (int jy = 0; jy < local_N_y; jy++) {
          int global_ix = local_start_x + ix;
          int global_jy = local_start_y + jy;

          field_data.index_i       = global_ix;
          field_data.index_j       = global_jy;
          field_data.index_k       = 0;
          field_data.beam_index    = di;
          field_data.N_frequencies = N_frequencies;
          field_data.stokes_I      = (THDF_float_t *)malloc(N_frequencies * sizeof(THDF_float_t));
          field_data.stokes_QI     = (THDF_float_t *)malloc(N_frequencies * sizeof(THDF_float_t));
          field_data.stokes_UI     = (THDF_float_t *)malloc(N_frequencies * sizeof(THDF_float_t));
          field_data.stokes_VI     = (THDF_float_t *)malloc(N_frequencies * sizeof(THDF_float_t));

          // Fill with example data
          for (int f = 0; f < N_frequencies; f++) {
            field_data.stokes_I[f] = sin(1.0 * (global_ix + 1) / (0.2 * N_x) * (global_jy + 1) * (f + 1) / N_frequencies);
            field_data.stokes_QI[f] =
                cos(1.5 * (global_ix + 1) / (0.2 * N_x) * (global_jy + 1) * (f + 1) / N_frequencies);
            field_data.stokes_UI[f] =
                sin(2.0 * (global_ix + 1) / (0.2 * N_x) * (global_jy + 1) * (f + 1) / N_frequencies);
            field_data.stokes_VI[f] =
                cos(3.0 * (global_ix + 1) / (0.2 * N_x) * (global_jy + 1) * (f + 1) / N_frequencies);
          }

          // Write field data for this position
          if (THDF_write_ard_field_to_hdf5(dset_handler, &field_data, global_ix, global_jy, 1, 1, N_frequencies) != 0) {
            fprintf(stderr, "Rank %d: Error writing ARD field data at (%d, %d) for direction %d\n", mpi_rank, global_ix,
                    global_jy, di);
            MPI_Abort(MPI_COMM_WORLD, -1);
          }

          free(field_data.stokes_I);
          free(field_data.stokes_QI);
          free(field_data.stokes_UI);
          free(field_data.stokes_VI);
        }
      }

      THDF_close_ard_field_handler_mpi(dset_handler);
    }

    H5Fclose(file_id);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_free(&writer_comm);
  MPI_Finalize();

  free(frequencies);
  free(muv);
  free(chiv);

  return 0;
}

/////////////////////////////////////////////////////////////
/// main
/////////////////////////////////////////////////////////////
int
main(int argc, char **argv) {
  // Switch between different examples:
  // return main_default(argc, argv);              // Simple MPI example
  // return main_domain_dec(argc, argv);  // Domain decomposition for regular fields
  // return main_example_ar_domain_dec(argc, argv);  // Domain decomposition for ARD fields
  return main_JKQ_domain_dec(argc, argv);  // Domain decomposition for JKQ fields
}