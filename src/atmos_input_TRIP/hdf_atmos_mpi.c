#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf_atmos.h"
#include "hdf_atmos_continuum.h"

// Compile with: mpicc -o hdf_atmos_mpi_example hdf_atmos_mpi.c
// hdf_atmos_functions.c -lhdf5 -lhdf5_hl

void
make_data_for_continuum(int N_frequencies, const int index_i, const int index_j, const int index_k,    //
                        double *sigma_c, double *epsilon_c, double *c_therm_emissivity_epsilon_c) {  //
  for (int i = 0; i < N_frequencies; i++) {
    sigma_c[i]                        = sin(0.1 * i + index_i + index_j + index_k) + 1.0;         //
    epsilon_c[i]                      = cos(0.05 * i + index_j + index_k) + 1.0;                  //
    c_therm_emissivity_epsilon_c[i] = sin(0.07 * i + index_k) + cos(0.03 * i + index_i) + 1.0;  //
  }
}

int
read_main_demo_cont() {

  int exit_code = EXIT_SUCCESS;

  char file_name[] = "atmosphere_mpi.h5";

  hid_t file = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    exit_code = EXIT_FAILURE;
    fprintf(stderr, "Error opening HDF5 file '%s'\n", file_name);
    goto finalize;
  }

  int     N_x, N_y, N_z;
  double *heights = NULL;
  if (THDF_read_geometry_from_hdf5(file, &N_x, &N_y, &N_z, &heights) < 0) {
    exit_code = EXIT_FAILURE;
    fprintf(stderr, "Error reading geometry data\n");
    goto free_heights;
  }

  printf("Read geometry: N_x=%d, N_y=%d, N_z=%d\n", N_x, N_y, N_z);
  printf("Heights:\n");
  for (int k = 0; k < N_z; k++) {
    printf("  height[%d] = %e cm, ", k, heights[k]);
  }
  printf("\n");

  int     N_frequencies;
  double *frequencies = (double *)malloc(N_frequencies * sizeof(double));
  if (THDF_read_frequency_grid_from_hdf5(file, &N_frequencies, &frequencies) < 0) {
    exit_code = EXIT_FAILURE;
    fprintf(stderr, "Error reading frequency grid data\n");
    goto free_frequencies;
  }

  printf("Read frequency grid: N_frequencies=%d\n", N_frequencies);
  printf("Frequencies:\n");
  for (int i = 0; i < N_frequencies; i++) {
    printf("  frequency[%d] = %e Hz, ", i, frequencies[i]);
  }
  printf("\n");

  THDF_continuum_t cont;
  double          *sigma_c                        = (double *)malloc(N_frequencies * sizeof(double));
  double          *epsilon_c                      = (double *)malloc(N_frequencies * sizeof(double));
  double          *c_therm_emissivity_epsilon_c = (double *)malloc(N_frequencies * sizeof(double));

  if (read_continuum_from_hdf5(file, N_frequencies, 4, 4, 4, 1, 1, 1, &cont, sigma_c, epsilon_c,
                               c_therm_emissivity_epsilon_c) < 0) {
    fprintf(stderr, "Error reading continuum data\n");
    exit_code = EXIT_FAILURE;
    goto free_sigma_c;
  }

  printf("Read continuum data for grid point (0,0,0):\n");
  for (int i = 0; i < N_frequencies; i++) {
    printf("  freq=%e Hz, sigma_c=%e, epsilon_c=%e, c_therm_emissivity_epsilon_c=%e\n", frequencies[i], sigma_c[i],
           epsilon_c[i], c_therm_emissivity_epsilon_c[i]);
  }

  THDF_atmos_t atmos_data;
  THDF_read_atmos_3D_point_hdf5(file, &atmos_data, 5, 5, 5);

  THDF_print_atmos_point(&atmos_data);

free_sigma_c:
  if (sigma_c) {
    free(sigma_c);
  }
  if (epsilon_c) {
    free(epsilon_c);
  }
  if (c_therm_emissivity_epsilon_c) {
    free(c_therm_emissivity_epsilon_c);
  }

free_heights:
  if (heights) {
    free(heights);
  }

free_frequencies:
  if (frequencies) {
    free(frequencies);
  }

finalize:
  H5Fclose(file);
  return exit_code;
}

#ifndef HDF_ATMOS_BUILD_LIBRARY
int
main(int argc, char **argv) {
  int mpi_rank, mpi_size;

  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  hid_t file, plist_id;

  // Default values
  int N_x           = 128;
  int N_y           = 128;
  int N_z           = 128;
  int N_frequencies = 96;

  // Parse command line arguments
  if (argc >= 5) {
    N_x           = atoi(argv[1]);
    N_y           = atoi(argv[2]);
    N_z           = atoi(argv[3]);
    N_frequencies = atoi(argv[4]);

    if (mpi_rank == 0) {
      printf("Using command line arguments: N_x=%d, N_y=%d, N_z=%d, N_frequencies=%d\n", N_x, N_y, N_z, N_frequencies);
    }
  } else if (mpi_rank == 0) {
    printf("Usage: %s <N_x> <N_y> <N_z> <N_frequencies>\n", argv[0]);
    printf("Using default values: N_x=%d, N_y=%d, N_z=%d, N_frequencies=%d\n", N_x, N_y, N_z, N_frequencies);
  }

  double *heights     = (double *)malloc(N_z * sizeof(double));
  double *frequencies = (double *)malloc(N_frequencies * sizeof(double));

  // Initialize example data
  for (int k = 0; k < N_z; k++) {
    heights[k] = k * 1.0e5;  // Example: height in cm
  }
  for (int i = 0; i < N_frequencies; i++) {
    frequencies[i] = 5.0e14 + i * 1.0e11;  // Example: frequency in Hz
  }

  // Set up file access property list with parallel MPI-IO access
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

  // Create the file collectively (all processes participate)
  file = H5Fcreate("atmosphere_mpi.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  if (file < 0) {
    if (mpi_rank == 0) {
      fprintf(stderr, "Error creating HDF5 file with MPI\n");
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  // Synchronize all processes
  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    printf("File created, all ranks ready\n");
    fflush(stdout);
  }

  // Close the parallel file and let rank 0 write metadata with serial I/O
  H5Fclose(file);
  MPI_Barrier(MPI_COMM_WORLD);

  // Example atom data (all processes need this)
  THDF_atom_two_levels_t atom = {
      .atomic_number = 1,
      .atomic_mass   = 1.0,
      .E_lower       = 10.0,
      .E_upper       = 20.0,
      .g_lower       = 2.0,
      .g_upper       = 4.0,
      .jl2           = 1,
      .ju2           = 2,
      .Aul           = 6.3e8,
      .a_coef_D2     = 1.0e-9,
      .b_coef_D2     = 0.4,
  };

  // Example two-term atom data (all processes need this)
  THDF_atom_two_terms_t atom_two_terms = {
      .atomic_number = 12,        // Example: Magnesium
      .atomic_mass   = 24.305,    // Atomic mass of Magnesium in atomic mass units
      .Aul           = 3.37e7,    // Einstein coefficient.
      .L2_upper      = 0,         // Upper term orbital angular momentum (multiplied by 2)
      .L2_lower      = 2,         // Lower term orbital angular momentum (multiplied by 2)
      .S2            = 2,         // Total spin (multiplied by 2)
      .E_size        = 3,         // Number of energy levels for this term
      .E_lower       = {21850.405, 21870.464, 21911.178},
      .E_upper       = {35051.264, 35051.264, 35051.264},
  };

  // Only rank 0 writes metadata using serial HDF5
  if (mpi_rank == 0) {
    file = H5Fopen("atmosphere_mpi.h5", H5F_ACC_RDWR, H5P_DEFAULT);

    if (write_geometry_to_hdf5(file, N_x, N_y, N_z, heights) < 0) {
      fprintf(stderr, "Error writing geometry data\n");
    }

    if (write_frequency_grid_to_hdf5(file, N_frequencies, frequencies) < 0) {
      fprintf(stderr, "Error writing frequency grid data\n");
    }

    if (write_atom_to_hdf5(file, &atom) < 0) {
      fprintf(stderr, "Error writing atom data\n");
    }

    if (write_atom_two_terms_to_hdf5(file, &atom_two_terms) < 0) {
      fprintf(stderr, "Error writing two-term atom data\n");
    }

    make_continuum_hdf5_mpi(file, N_x, N_y, N_z, N_frequencies);  // Example dimensions
    THDF_create_continuum_4D_matrices(file, N_x, N_y, N_z, N_frequencies);

    H5Fclose(file);
    printf("Rank 0: Metadata written successfully\n");
    fflush(stdout);
  }  // END MPI rank 0 metadata writing

  // Synchronize all processes
  MPI_Barrier(MPI_COMM_WORLD);
  double start_time = MPI_Wtime();

  if (mpi_rank == 0) {
    printf("Rank 0: Creating atmosphere dataset collectively...\n");
    fflush(stdout);
  }

  // Calculate work distribution for each MPI process
  int local_N_x = N_x / mpi_size;
  int remainder = N_x % mpi_size;
  int start_i   = mpi_rank * local_N_x + (mpi_rank < remainder ? mpi_rank : remainder);
  if (mpi_rank < remainder)
    local_N_x++;

  // Allocate local atmosphere data
  THDF_atmos_t *local_atmos_data = (THDF_atmos_t *)malloc(local_N_x * N_y * N_z * sizeof(THDF_atmos_t));
  if (local_atmos_data == NULL) {
    fprintf(stderr, "Rank %d: Error allocating memory for local atmosphere data\n", mpi_rank);
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  // All processes reopen file with parallel I/O
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file = H5Fopen("atmosphere_mpi.h5", H5F_ACC_RDWR, plist_id);
  H5Pclose(plist_id);

  if (file < 0) {
    fprintf(stderr, "Rank %d: Error reopening file for parallel I/O\n", mpi_rank);
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  // Create MPI-enabled atmosphere dataset (collective operation)
  THDF_atmos_dataset_handler_t *atmos_dset = create_atmosphere_dataset_mpi(file, N_x, N_y, N_z, H5P_DEFAULT);

  if (atmos_dset == NULL) {
    fprintf(stderr, "Rank %d: Error creating MPI atmosphere dataset\n", mpi_rank);
    free(local_atmos_data);
    H5Fclose(file);
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  // Placeholder for the normalized output datasets.
  // In this case local_N_x, N_y, N_z is the local size of the slice
  // howned by the current MPI process.
  // in this example the domain in subdivided of vertical slices along the x-axis.
  u_int16_t *norm_sigma_c   = (uint16_t *)malloc(N_frequencies * local_N_x * N_y * N_z * sizeof(uint16_t));
  u_int16_t *norm_epsilon_c = (uint16_t *)malloc(N_frequencies * local_N_x * N_y * N_z * sizeof(uint16_t));
  u_int16_t *norm_c_therm_emissivity_epsilon_c =
      (uint16_t *)malloc(N_frequencies * local_N_x * N_y * N_z * sizeof(uint16_t));

  // Input continuum dataset.
  double           *sigma_c                        = (double *)malloc(N_frequencies * sizeof(double));
  double           *epsilon_c                      = (double *)malloc(N_frequencies * sizeof(double));
  double           *c_therm_emissivity_epsilon_c = (double *)malloc(N_frequencies * sizeof(double));
  THDF_continuum_t *cont = (THDF_continuum_t *)malloc(local_N_x * N_y * N_z * sizeof(THDF_continuum_t));

  // Fill local data with MPI-aware atmospheric data
  for (int i = 0; i < local_N_x; i++) {
    for (int j = 0; j < N_y; j++) {
      for (int k = 0; k < N_z; k++) {

        // Generate a fake data set and writes it in sigma_c, epsilon_c, c_therm_emissivity_epsilon_c
        make_data_for_continuum(N_frequencies, start_i + i, j, k, sigma_c, epsilon_c, c_therm_emissivity_epsilon_c);

        const int offset_cont = (i * N_y * N_z + j * N_z + k) * N_frequencies;

        int local_idx   = i * N_y * N_z + j * N_z + k;
        cont[local_idx] = normalize_continuum_data(sigma_c,                                             //
                                                   epsilon_c,                                           //
                                                   c_therm_emissivity_epsilon_c,                      //
                                                   N_frequencies,                                       //
                                                   &norm_sigma_c[offset_cont],                          //
                                                   &norm_epsilon_c[offset_cont],                        //
                                                   &norm_c_therm_emissivity_epsilon_c[offset_cont]);  //

        // Set the current absolute index in the cont-struct
        int global_i = start_i + i;

        cont[local_idx].index_i = global_i;
        cont[local_idx].index_j = j;
        cont[local_idx].index_k = k;

        // genertae a fake THDF_atmos_t data for this point.
        // And copy it in the output placeholder.
        local_atmos_data[local_idx] = create_atmos_data_point_mpi(global_i, j, k, mpi_rank, &atom);
      }
    }
  }

  // Set up collective write property list
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,
                   H5FD_MPIO_INDEPENDENT);  // Use independent mode for now

  // Write local data collectively
  int result = write_atmosphere_data_to_dataset_mpi(atmos_dset,                      //
                                                    local_atmos_data,                //
                                                    start_i, 0, 0,                   //
                                                    local_N_x, N_y, N_z, plist_id);  //

  // Write the continuum
  int result_cont =
      write_normalized_continuum_to_hdf5(file, N_frequencies,                                                       //
                                         start_i, 0, 0,                                                             //
                                         local_N_x, N_y, N_z,                                                       //
                                         cont, norm_sigma_c, norm_epsilon_c, norm_c_therm_emissivity_epsilon_c);  //

  if (result < 0) {
    fprintf(stderr, "Rank %d: Error writing local atmosphere data\n", mpi_rank);
  }

  if (result_cont < 0) {
    fprintf(stderr, "Rank %d: Error writing local continuum data\n", mpi_rank);
  }

  // Clean up
  close_atmosphere_dataset_mpi(atmos_dset);
  free(local_atmos_data);

  H5Pclose(plist_id);
  H5Fclose(file);

  double end_time = MPI_Wtime();

  if (mpi_rank == 0) {
    printf("MPI HDF5 atmosphere file created successfully with %d processes\n\n", mpi_size);

    // Rank 0 verifies the file by reading it back
    printf("Rank 0: Verifying the MPI-generated file...\n\n");
    int verify_result = THDF_read_main_demo_from_file("atmosphere_mpi.h5", false);

    if (verify_result == EXIT_SUCCESS) {
      printf("\nRank 0: File verification successful!\n");
    } else {
      printf("\nRank 0: File verification failed!\n");
    }

    printf("Time taken for parallel write: %e seconds\n", end_time - start_time);
  }

  free(heights);
  free(frequencies);
  free(sigma_c);
  free(epsilon_c);
  free(c_therm_emissivity_epsilon_c);
  free(norm_sigma_c);
  free(norm_epsilon_c);
  free(norm_c_therm_emissivity_epsilon_c);
  free(cont);

  if (mpi_rank == 0) {
    // How to read ....
    read_main_demo_cont();
  }

  MPI_Finalize();

  return EXIT_SUCCESS;
}
#endif /* HDF_ATMOS_BUILD_LIBRARY */