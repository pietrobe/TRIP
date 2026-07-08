#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf_atmos.h"

//============================================================================
// read_main_demo_from_file_TRIP
//============================================================================

int
THDF_read_main_demo_from_file_TRIP(const char *filename, const bool verbose) {

  hid_t file;

  printf("==== Reading HDF5 file: %s ====\n\n", filename);

  file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    fprintf(stderr, "Error opening HDF5 file\n");
    return EXIT_FAILURE;
  }

  int     N_x, N_y, N_z;
  double *heights = NULL;
  if (THDF_read_geometry_from_hdf5(file, &N_x, &N_y, &N_z, &heights) < 0) {
    fprintf(stderr, "Error reading geometry data\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }

  printf("Geometry: N_x=%d, N_y=%d, N_z=%d\n", N_x, N_y, N_z);
  printf("Heights:\n");

  if (verbose) {

    for (int k = 0; k < N_z; k++) {
      printf("  height[%d] = %e cm\n", k, heights[k]);
    }
  }
  free(heights);

  THDF_atom_two_levels_t atom;
  if (THDF_read_atom_from_hdf5(file, &atom) != EXIT_SUCCESS) {
    fprintf(stderr, "Error reading atom data\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }

  printf("Atom data:\n");
  printf("  atomic_mass = %f\n", atom.atomic_mass);
  printf("  E_lower = %f\n", atom.E_lower);
  printf("  E_upper = %f\n", atom.E_upper);
  printf("  g_lower = %f\n", atom.g_lower);
  printf("  g_upper = %f\n", atom.g_upper);
  printf("  jl2 = %d\n", atom.jl2);
  printf("  ju2 = %d\n", atom.ju2);
  printf("  Aul = %e\n", atom.Aul);

  THDF_atmos_t *atmos_data = NULL;
  if (THDF_read_atmos_3D_data_from_hdf5(file, &atmos_data, &N_x, &N_y, &N_z) != EXIT_SUCCESS) {
    fprintf(stderr, "Error reading atmosphere data\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }

  THDF_read_atmos_3D_data_from_hdf5(file, &atmos_data, &N_x, &N_y, &N_z);

  for (int k = 0; k < N_z; k++) {
    for (int j = 0; j < N_y; j++) {
      for (int i = 0; i < N_x; i++) {
        int          idx   = i * N_y * N_z + j * N_z + k;
        THDF_atmos_t point = atmos_data[idx];
        if (point.index_i != i || point.index_j != j || point.index_k != k) {
          fprintf(stderr, "Data index mismatch at (%d, %d, %d)\n", i, j, k);
          free(atmos_data);
          H5Fclose(file);
          return EXIT_FAILURE;
        }
      }
    }
  }

  printf("Atmosphere data point at (i=0, j=0, k=0):\n");
  THDF_print_atmos_point(&atmos_data[33]);

  free(atmos_data);
  H5Fclose(file);

  return EXIT_SUCCESS;
}

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
write_main_demo_TRIP(int argc, char **argv) {

  (void)argc;
  (void)argv;

  // HDF5 file identifier
  hid_t file;

  // Define the spatial grid dimensions.
  int N_x = 30;
  int N_y = 50;
  int N_z = 20;

  file = H5Fcreate("atmosphere.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    fprintf(stderr, "Error creating HDF5 file\n");
    return EXIT_FAILURE;
  }

  // Write geometry data to the HDF5 file

  double *heights = (double *)malloc(N_z * sizeof(double));
  if (heights == NULL) {
    fprintf(stderr, "Error allocating memory for geometry heights\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }
  for (int k = 0; k < N_z; k++) {
    heights[k] = k * 1.0e5;  // Example: height in cm
  }

  int geometry_flag = write_geometry_to_hdf5(file, N_x, N_y, N_z, heights);

  if (geometry_flag < 0) {
    fprintf(stderr, "Error writing geometry data\n");
    free(heights);
    H5Fclose(file);
    return EXIT_FAILURE;
  }
  free(heights);

  // Example atom data
  THDF_atom_two_levels_t atom;
  // Example for Ca I 4227 Å line
  atom.atomic_number = 20;         // Example: Calcium
  atom.atomic_mass   = 40.078;     // Atomic mass of Calcium in atomic mass units
  atom.E_lower       = 0.0;        // Ground state energy in erg
  atom.E_upper       = 23652.304;  // Excited state energy in erg
  atom.g_lower       = 0.0;        // g lower level Lande factors
  atom.g_upper       = 1.0;        // g upper level Lande factors
  atom.jl2           = 0;          // Total angular momentum of lower level
  atom.ju2           = 2;          // Total angular momentum of upper level
  atom.Aul           = 2.18E+08;   // Einstein coefficient.

  // Write atom data to HDF5 file
  int flag_atom = write_atom_to_hdf5(file, &atom);

  // Check for errors
  if (flag_atom < 0) {
    fprintf(stderr, "Error writing atom data\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }

  // Example two-term atom data
  THDF_atom_two_terms_t atom_two_terms;
  // Example for Mg I b2 triplet (3s3p 3P - 3s4s 3S transition)
  atom_two_terms.atomic_number = 12;        // Example: Magnesium
  atom_two_terms.atomic_mass   = 24.305;    // Atomic mass of Magnesium in atomic mass units
  atom_two_terms.Aul           = 3.37E+07;  // Einstein coefficient.
  atom_two_terms.L2_upper      = 0;         // Upper term orbital angular momentum (multiplied by 2)
  atom_two_terms.L2_lower      = 2;         // Lower term orbital angular momentum (multiplied by 2)
  atom_two_terms.S2            = 2;         // Total spin (multiplied by 2)

  atom_two_terms.E_size     = 3;  // Number of energy levels for this term
  atom_two_terms.E_lower[0] = 21850.405;
  atom_two_terms.E_lower[1] = 21870.464;
  atom_two_terms.E_lower[2] = 21911.178;
  atom_two_terms.E_upper[0] = 35051.264;
  atom_two_terms.E_upper[1] = 35051.264;
  atom_two_terms.E_upper[2] = 35051.264;

  // Write two-term atom data to HDF5 file
  int flag_atom_two_terms = write_atom_two_terms_to_hdf5(file, &atom_two_terms);

  // Check for errors
  if (flag_atom_two_terms < 0) {
    fprintf(stderr, "Error writing two-term atom data\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }

  // Allocate a sufficiently large array to hold the atmosphere data for the actual exigence.
  THDF_atmos_t *atmos_data = (THDF_atmos_t *)malloc(N_x * N_y * N_z * sizeof(THDF_atmos_t));
  if (atmos_data == NULL) {
    fprintf(stderr, "Error allocating memory for atmosphere data\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }

  // Open and generte the entry in the HDF5 file for the atmosphere dataset.
  THDF_atmos_dataset_handler_t *atmos_dset = THDF_create_atmosphere_dataset(file, N_x, N_y, N_z);

  // It is important to follow this loop order (i, j, k) if you want to write the data in blocks.
  // In this case we write one block of size (1, 1, N_z) for each (i, j) pair.
  // But if you want you can write bloks for any sizes from (1, 1, 1) to (N_x, N_y, N_z)
  // with the appropriate loop order and block sizes.
  // IMPORTNT: This method wors only in serial. Parallel (MPI, OpenMP, Thareds) writing requires a different approach.
  for (int i = 0; i < N_x; i++) {
    for (int j = 0; j < N_y; j++) {

      int idx_base = i * N_y * N_z + j * N_z;

      for (int k = 0; k < N_z; k++) {

        int idx         = idx_base + k;
        atmos_data[idx] = create_atmos_data_point(i, j, k, &atom);
      }
      THDF_write_atmosphere_data_to_dataset(atmos_dset,             //
                                            &atmos_data[idx_base],  // Pointer to the first point of the block (i, j, 0)
                                            i, j, 0,                // Start indices
                                            1, 1, N_z);             // Count (block) sizes
    }
  }

  // Close the atmosphere dataset after writing all the data.
  THDF_close_atmosphere_dataset(atmos_dset);
  H5Fclose(file);

  //   if (write_atmos_3D_data_to_hdf5(file, atmos_data, N_x, N_y, N_z) !=
  //       EXIT_SUCCESS) {
  //     fprintf(stderr, "Error writing atmosphere data\n");
  //     free(atmos_data);
  //     H5Fclose(file);
  //     return EXIT_FAILURE;
  //   }

  free(atmos_data);
  atmos_data = NULL;

  printf("\nFile 'atmosphere.h5' created successfully!\n\n");

  return EXIT_SUCCESS;
}

int
main(int argc, char **argv) {

  // Check if launched with multiple MPI ranks (e.g. via mpirun)
  // This is a serial example, so we want to prevent accidental parallel launches
  // If your application is serial, you should use this as example.
  int launcher_size = get_launcher_world_size();
  if (launcher_size > 1) {
    fprintf(stderr,
            "This is the serial executable and should not be launched with multiple MPI ranks.\n"
            "Detected launcher world size = %d.\n",
            launcher_size);
    fprintf(stderr, "Use './hdf_atmos_example' for serial runs or 'mpirun ./hdf_atmos_mpi_example' for parallel runs.\n");
    return EXIT_FAILURE;
  }

  // First, create and write the file
  write_main_demo_TRIP(argc, argv);

  // Then read it back
  return read_main_demo();
}
