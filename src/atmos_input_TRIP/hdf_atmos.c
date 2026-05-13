#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf_atmos.h"
#include "hdf_atmos_continuum.h"

// Compile with: gcc -o hdf_atmos_example hdf_atmos.c -lhdf5 -lhdf5_hl

// Forward declarations
int
read_main_demo();
int
write_main_demo(int argc, char **argv);

static int
get_launcher_world_size(void) {
  const char *size_vars[] = {"OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMIX_SIZE", "MV2_COMM_WORLD_SIZE"};
  const int   n_vars      = (int)(sizeof(size_vars) / sizeof(size_vars[0]));

  for (int idx = 0; idx < n_vars; idx++) {
    const char *value = getenv(size_vars[idx]);
    if (value != NULL && value[0] != '\0') {
      int parsed_size = atoi(value);
      if (parsed_size > 0) {
        return parsed_size;
      }
    }
  }

  return 1;
}

////////////////////////////////////////////////////////////////////

int
read_main_demo() {
  return THDF_read_main_demo_from_file("atmosphere.h5", true);
}

int
write_main_demo(int argc, char **argv) {

  (void)argc;
  (void)argv;

  hid_t file;

  int N_x = 30;
  int N_y = 50;
  int N_z = 20;

  double heights[20];
  double frequencies[50];
  for (int k = 0; k < N_z; k++) {
    heights[k] = k * 1.0e5;  // Example: height in cm
  }
  for (int i = 0; i < 50; i++) {
    frequencies[i] = 5.0e14 + i * 1.0e11;  // Example: frequency in Hz
  }

  file = H5Fcreate("atmosphere.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    fprintf(stderr, "Error creating HDF5 file\n");
    return EXIT_FAILURE;
  }

  // Write geometry data to the HDF5 file
  if (write_geometry_to_hdf5(file, N_x, N_y, N_z, heights) < 0) {
    fprintf(stderr, "Error writing geometry data\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }

  // Write frequency grid data to the HDF5 file
  if (write_frequency_grid_to_hdf5(file, 50, frequencies) < 0) {
    fprintf(stderr, "Error writing frequency grid data\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }

  // Example atom data
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
  if (write_atom_to_hdf5(file, &atom) < 0) {
    fprintf(stderr, "Error writing atom data\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }

  THDF_atmos_t *atmos_data = (THDF_atmos_t *)malloc(N_x * N_y * N_z * sizeof(THDF_atmos_t));
  if (atmos_data == NULL) {
    fprintf(stderr, "Error allocating memory for atmosphere data\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }

  THDF_atmos_dataset_handler_t *atmos_dset = THDF_create_atmosphere_dataset(file, N_x, N_y, N_z);

  for (int i = 0; i < N_x; i++) {
    for (int j = 0; j < N_y; j++) {

      int idx_base = i * N_y * N_z + j * N_z;

      for (int k = 0; k < N_z; k++) {

        int idx         = idx_base + k;
        atmos_data[idx] = create_atmos_data_point(i, j, k, &atom);
      }
      THDF_write_atmosphere_data_to_dataset(atmos_dset, &atmos_data[idx_base], i, j, 0, 1, 1, N_z);
    }
  }

  THDF_close_atmosphere_dataset(atmos_dset);

  //   if (write_atmos_3D_data_to_hdf5(file, atmos_data, N_x, N_y, N_z) !=
  //       EXIT_SUCCESS) {
  //     fprintf(stderr, "Error writing atmosphere data\n");
  //     free(atmos_data);
  //     H5Fclose(file);
  //     return EXIT_FAILURE;
  //   }

  free(atmos_data);
  atmos_data = NULL;

  H5Fclose(file);

  printf("\nFile 'atmosphere.h5' created successfully!\n\n");

  return EXIT_SUCCESS;
}

int
main(int argc, char **argv) {

  return continuum_main_demo(argc, argv);

  int launcher_size = get_launcher_world_size();
  if (launcher_size > 1) {
    fprintf(stderr,
            "This is the serial executable and should not be launched with multiple MPI ranks.\n"
            "Detected launcher world size = %d.\n"
            "Use './hdf_atmos_example' for serial runs or 'mpirun ./hdf_atmos_mpi_example' for parallel runs.\n",
            launcher_size);
    return EXIT_FAILURE;
  }

  // First, create and write the file
  write_main_demo(argc, argv);

  // Then read it back
  return read_main_demo();
}
