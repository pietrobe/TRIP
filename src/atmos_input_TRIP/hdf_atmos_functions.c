#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf_atmos.h"

// Compile with: gcc -o hdf_atmos_example hdf_atmos.c -lhdf5 -lhdf5_hl

//============================================================================
// create_hdf_atmos_compound_type
//============================================================================

hid_t
create_hdf_atmos_compound_type(void) {
  hid_t atmos_type = H5Tcreate(H5T_COMPOUND, sizeof(THDF_atmos_t));

  H5Tinsert(atmos_type, "index_i", HOFFSET(THDF_atmos_t, index_i), H5T_NATIVE_INT16);
  H5Tinsert(atmos_type, "index_j", HOFFSET(THDF_atmos_t, index_j), H5T_NATIVE_INT16);
  H5Tinsert(atmos_type, "index_k", HOFFSET(THDF_atmos_t, index_k), H5T_NATIVE_INT16);
  H5Tinsert(atmos_type, "temperature", HOFFSET(THDF_atmos_t, temperature), H5T_NATIVE_FLOAT);
  H5Tinsert(atmos_type, "vmicro", HOFFSET(THDF_atmos_t, vmicro), H5T_NATIVE_DOUBLE);
  H5Tinsert(atmos_type, "Doppler_width", HOFFSET(THDF_atmos_t, Doppler_width), H5T_NATIVE_DOUBLE);
  H5Tinsert(atmos_type, "damping", HOFFSET(THDF_atmos_t, damping), H5T_NATIVE_DOUBLE);
  H5Tinsert(atmos_type, "pop_lower_level", HOFFSET(THDF_atmos_t, pop_lower_level), H5T_NATIVE_DOUBLE);
  H5Tinsert(atmos_type, "pop_upper_level", HOFFSET(THDF_atmos_t, pop_upper_level), H5T_NATIVE_DOUBLE);
  H5Tinsert(atmos_type, "Cul", HOFFSET(THDF_atmos_t, Cul), H5T_NATIVE_DOUBLE);
  H5Tinsert(atmos_type, "Qel", HOFFSET(THDF_atmos_t, Qel), H5T_NATIVE_DOUBLE);
  H5Tinsert(atmos_type, "nH", HOFFSET(THDF_atmos_t, nH), H5T_NATIVE_DOUBLE);
  H5Tinsert(atmos_type, "D2", HOFFSET(THDF_atmos_t, D2), H5T_NATIVE_DOUBLE);
  H5Tinsert(atmos_type, "magnetic_field_x", HOFFSET(THDF_atmos_t, magnetic_field_x), H5T_NATIVE_FLOAT);
  H5Tinsert(atmos_type, "magnetic_field_y", HOFFSET(THDF_atmos_t, magnetic_field_y), H5T_NATIVE_FLOAT);
  H5Tinsert(atmos_type, "magnetic_field_z", HOFFSET(THDF_atmos_t, magnetic_field_z), H5T_NATIVE_FLOAT);
  H5Tinsert(atmos_type, "bulk_velocity_x", HOFFSET(THDF_atmos_t, bulk_velocity_x), H5T_NATIVE_FLOAT);
  H5Tinsert(atmos_type, "bulk_velocity_y", HOFFSET(THDF_atmos_t, bulk_velocity_y), H5T_NATIVE_FLOAT);
  H5Tinsert(atmos_type, "bulk_velocity_z", HOFFSET(THDF_atmos_t, bulk_velocity_z), H5T_NATIVE_FLOAT);
  H5Tinsert(atmos_type, "c_scat_opacity_sigma_c", HOFFSET(THDF_atmos_t, c_scat_opacity_sigma_c), H5T_NATIVE_DOUBLE);
  H5Tinsert(atmos_type, "c_therm_opacity_k_c", HOFFSET(THDF_atmos_t, c_therm_opacity_k_c), H5T_NATIVE_DOUBLE);
  H5Tinsert(atmos_type, "c_therm_emissivity_epsilon_c", HOFFSET(THDF_atmos_t, c_therm_emissivity_epsilon_c),
            H5T_NATIVE_DOUBLE);

  return atmos_type;
}

//============================================================================
// create_atmos_data_point
//============================================================================

THDF_atmos_t
create_atmos_data_point(int i, int j, int k, const THDF_atom_two_levels_t *atom) {
  THDF_atmos_t atmos_point;

  atmos_point.index_i                      = i;
  atmos_point.index_j                      = j;
  atmos_point.index_k                      = k;
  atmos_point.temperature                  = 5000.0 + k * 100.0 + i * 10.0 + j * 5.0;
  atmos_point.vmicro                       = 2.0 + 0.1 * k;
  atmos_point.Doppler_width                = 1.0e8 + 1.0e6 * k;
  atmos_point.damping                      = 1.0e-4 + 1.0e-6 * k;
  atmos_point.pop_lower_level              = 1.0e12 - k * 1.0e10;
  atmos_point.pop_upper_level              = 1.0e10 + k * 1.0e9;
  atmos_point.Cul                          = 1.0e-8 + 1.0e-10 * k;
  atmos_point.Qel                          = 1.0e-9 + 1.0e-11 * k;
  atmos_point.nH                           = 1.0e17 + 1.0e15 * k;
  atmos_point.magnetic_field_x             = (float)(100.0 + 10.0 * k);
  atmos_point.magnetic_field_y             = (float)(45.0 + k);
  atmos_point.magnetic_field_z             = (float)(90.0 + 2 * k);
  atmos_point.bulk_velocity_x              = (float)(1.0e5 + 1.0e4 * k);
  atmos_point.bulk_velocity_y              = (float)(0.0 + k);
  atmos_point.bulk_velocity_z              = (float)(0.0 + 2 * k);
  atmos_point.c_scat_opacity_sigma_c       = 1.0e-2 + 1.0e-4 * k;
  atmos_point.c_therm_opacity_k_c          = 1.0e-3 + 1.0e-5 * k;
  atmos_point.c_therm_emissivity_epsilon_c = 1.0e-1 + 1.0e-3 * k;

  if (atom != NULL) {
    atmos_point.D2 = atom->a_coef_D2 * atmos_point.nH * pow(atmos_point.temperature / 5000.0, atom->b_coef_D2);
  } else {
    atmos_point.D2 = 0.0;
  }

  return atmos_point;
}

void
THDF_print_atmos_point(const THDF_atmos_t *point) {
  if (point == NULL) {
    printf("THDF_atmos_t: (null)\n");
    return;
  }

  printf("THDF_atmos_t {\n");
  printf("  index_i = %d\n", point->index_i);
  printf("  index_j = %d\n", point->index_j);
  printf("  index_k = %d\n", point->index_k);
  printf("  temperature = %f\n", point->temperature);
  printf("  vmicro = %e\n", point->vmicro);
  printf("  Doppler_width = %e\n", point->Doppler_width);
  printf("  damping = %e\n", point->damping);
  printf("  pop_lower_level = %e\n", point->pop_lower_level);
  printf("  pop_upper_level = %e\n", point->pop_upper_level);
  printf("  Cul = %e\n", point->Cul);
  printf("  Qel = %e\n", point->Qel);
  printf("  nH = %e\n", point->nH);
  printf("  magnetic_field_x = %f\n", point->magnetic_field_x);
  printf("  magnetic_field_y = %f\n", point->magnetic_field_y);
  printf("  magnetic_field_z = %f\n", point->magnetic_field_z);
  printf("  bulk_velocity_x = %f\n", point->bulk_velocity_x);
  printf("  bulk_velocity_y = %f\n", point->bulk_velocity_y);
  printf("  bulk_velocity_z = %f\n", point->bulk_velocity_z);
  printf("  c_scat_opacity_sigma_c = %e\n", point->c_scat_opacity_sigma_c);
  printf("  c_therm_opacity_k_c = %e\n", point->c_therm_opacity_k_c);
  printf("  c_therm_emissivity_epsilon_c = %e\n", point->c_therm_emissivity_epsilon_c);
  printf("  D2 = %e\n", point->D2);
  printf("}\n");
}

//============================================================================
// MPI-enabled functions
//============================================================================

//============================================================================
// create_atmosphere_dataset_mpi
//============================================================================
THDF_atmos_dataset_handler_t *
create_atmosphere_dataset_mpi(hid_t file, int N_x, int N_y, int N_z, hid_t plist_id) {
  THDF_atmos_dataset_handler_t *atmos_dset = (THDF_atmos_dataset_handler_t *)malloc(sizeof(THDF_atmos_dataset_handler_t));
  if (!atmos_dset) {
    fprintf(stderr, "Error allocating memory for atmosphere dataset struct\n");
    return NULL;
  }

  (void)plist_id;  // Currently unused, but can be used for collective dataset creation if needed

  // Initialize struct
  atmos_dset->file_id = file;
  atmos_dset->N_x     = N_x;
  atmos_dset->N_y     = N_y;
  atmos_dset->N_z     = N_z;
  atmos_dset->is_open = 0;

  hsize_t dims[3] = {N_x, N_y, N_z};

  atmos_dset->group_id = -1;  // dataset lives at root, no group

  // Create the compound datatype for hdf_atmos
  atmos_dset->datatype_id = create_hdf_atmos_compound_type();

  // Create dataspace and dataset with MPI support (all processes collectively)
  atmos_dset->dataspace_id = H5Screate_simple(3, dims, NULL);
  atmos_dset->dataset_id =
      H5Dcreate2(file, "atmos", atmos_dset->datatype_id, atmos_dset->dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (atmos_dset->dataset_id < 0) {
    fprintf(stderr, "Error creating atmosphere dataset with MPI\n");
    H5Tclose(atmos_dset->datatype_id);
    H5Sclose(atmos_dset->dataspace_id);
    free(atmos_dset);
    return NULL;
  }

  atmos_dset->is_open = 1;
  return atmos_dset;
}

//============================================================================
// write_atmosphere_data_to_dataset_mpi
//============================================================================
int
write_atmosphere_data_to_dataset_mpi(THDF_atmos_dataset_handler_t *atmos_dset, THDF_atmos_t *atmos_data, hsize_t start_i,
                                     hsize_t start_j, hsize_t start_k, hsize_t count_i, hsize_t count_j, hsize_t count_k,
                                     hid_t plist_id) {
  if (!atmos_dset || !atmos_dset->is_open) {
    fprintf(stderr, "Error: atmosphere dataset is not open\n");
    return -1;
  }

  // Define hyperslab for partial writing
  hsize_t start[3] = {start_i, start_j, start_k};
  hsize_t count[3] = {count_i, count_j, count_k};

  // Create memory space for the data to be written
  hid_t memspace = H5Screate_simple(3, count, NULL);

  // Select hyperslab in file dataspace
  H5Sselect_hyperslab(atmos_dset->dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

  // Write data collectively with MPI
  herr_t status =
      H5Dwrite(atmos_dset->dataset_id, atmos_dset->datatype_id, memspace, atmos_dset->dataspace_id, plist_id, atmos_data);

  H5Sclose(memspace);

  if (status < 0) {
    fprintf(stderr, "Error writing atmosphere data to dataset with MPI\n");
    return -1;
  }

  return 0;
}

//============================================================================
// close_atmosphere_dataset_mpi
//============================================================================
int
close_atmosphere_dataset_mpi(THDF_atmos_dataset_handler_t *atmos_dset) {
  if (!atmos_dset) {
    return -1;
  }

  if (atmos_dset->is_open) {
    H5Dclose(atmos_dset->dataset_id);
    H5Sclose(atmos_dset->dataspace_id);
    H5Tclose(atmos_dset->datatype_id);
    if (atmos_dset->group_id >= 0)
      H5Gclose(atmos_dset->group_id);
    atmos_dset->is_open = 0;
  }

  free(atmos_dset);
  return 0;
}

//============================================================================
// create_hdf_atmos_compound_type_mpi
//============================================================================
hid_t
create_hdf_atmos_compound_type_mpi(void) {
  // MPI version is identical to non-MPI version for compound types
  return create_hdf_atmos_compound_type();
}

//============================================================================
// create_atmos_data_point_mpi
//============================================================================
THDF_atmos_t
create_atmos_data_point_mpi(int i, int j, int k, int mpi_rank, const THDF_atom_two_levels_t *atom) {
  THDF_atmos_t atmos_point;

  atmos_point.index_i = i;
  atmos_point.index_j = j;
  atmos_point.index_k = k;

  // Add MPI rank dependency to some calculations for variety across processes
  atmos_point.temperature                  = 5000.0 + k * 100.0 + i * 10.0 + j * 5.0 + mpi_rank * 2.0;
  atmos_point.vmicro                       = 2.0 + 0.1 * k + mpi_rank * 0.01;
  atmos_point.Doppler_width                = 1.0e8 + 1.0e6 * k + mpi_rank * 1.0e5;
  atmos_point.damping                      = 1.0e-4 + 1.0e-6 * k + mpi_rank * 1.0e-7;
  atmos_point.pop_lower_level              = 1.0e12 - k * 1.0e10 + mpi_rank * 1.0e9;
  atmos_point.pop_upper_level              = 1.0e10 + k * 1.0e9 + mpi_rank * 1.0e8;
  atmos_point.Cul                          = 1.0 + 2.0 * k + mpi_rank;
  atmos_point.Qel                          = 1.0 + 3.0 * k + mpi_rank;
  atmos_point.nH                           = 1.0e17 + 1.0e15 * k + mpi_rank * 1.0e14;
  atmos_point.magnetic_field_x             = (float)(100.0 + 10.0 * k + mpi_rank * 1.0);
  atmos_point.magnetic_field_y             = (float)(45.0 + k + mpi_rank * 0.5);
  atmos_point.magnetic_field_z             = (float)(90.0 + 2 * k + mpi_rank * 1.0);
  atmos_point.bulk_velocity_x              = (float)(1.0e5 + 1.0e4 * k + mpi_rank * 1.0e3);
  atmos_point.bulk_velocity_y              = (float)(0.0 + k + mpi_rank * 0.1);
  atmos_point.bulk_velocity_z              = (float)(0.0 + 2 * k + mpi_rank * 0.2);
  atmos_point.c_scat_opacity_sigma_c       = 1.0e-2 + 1.0e-4 * k + mpi_rank * 1.0e-5;
  atmos_point.c_therm_opacity_k_c          = 1.0e-3 + 1.0e-5 * k + mpi_rank * 1.0e-6;
  atmos_point.c_therm_emissivity_epsilon_c = 1.0e-1 + 1.0e-3 * k + mpi_rank * 1.0e-4;

  if (atom != NULL) {
    atmos_point.D2 = atom->a_coef_D2 * atmos_point.nH * pow(atmos_point.temperature / 5000.0, atom->b_coef_D2);
  } else {
    atmos_point.D2 = 0.0;
  }

  return atmos_point;
}

//============================================================================
// write_geometry_to_hdf5
//============================================================================
int
write_geometry_to_hdf5(hid_t file, int N_x, int N_y, int N_z, double *heights) {

  THDF_atmos_geometry_t geometry = {N_z, N_x, N_y, NULL, 100000.0};

  geometry.height = (double *)malloc(N_z * sizeof(double));
  if (geometry.height == NULL) {
    fprintf(stderr, "Error allocating memory for height array\n");
    return -1;
  }

  for (int k = 0; k < N_z; k++) {
    geometry.height[k] = heights[k];
  }

  // Write geometry data to HDF5 file
  hid_t geom_group = H5Gcreate2(file, "/geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (geom_group < 0) {
    fprintf(stderr, "Error creating geometry group\n");
    free(geometry.height);
    return -1;
  }

  // Write N_height as scalar int64.
  hid_t scalar_space = H5Screate(H5S_SCALAR);
  if (scalar_space < 0) {
    fprintf(stderr, "Error creating scalar dataspace for geometry metadata\n");
    H5Gclose(geom_group);
    free(geometry.height);
    return -1;
  }

  int64_t n_height_i64 = (int64_t)geometry.N_z;
  hid_t   dset_n_height =
      H5Dcreate2(geom_group, "N_height", H5T_STD_I64LE, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset_n_height < 0 || H5Dwrite(dset_n_height, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n_height_i64) < 0) {
    fprintf(stderr, "Error writing geometry dataset N_height\n");
    if (dset_n_height >= 0)
      H5Dclose(dset_n_height);
    H5Sclose(scalar_space);
    H5Gclose(geom_group);
    free(geometry.height);
    return -1;
  }
  H5Dclose(dset_n_height);

  // Write N_ij as scalar int64.
  int64_t n_ij_i64  = (int64_t)geometry.N_x;
  hid_t   dset_n_ij = H5Dcreate2(geom_group, "N_ij", H5T_STD_I64LE, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset_n_ij < 0 || H5Dwrite(dset_n_ij, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n_ij_i64) < 0) {
    fprintf(stderr, "Error writing geometry dataset N_ij\n");
    if (dset_n_ij >= 0)
      H5Dclose(dset_n_ij);
    H5Sclose(scalar_space);
    H5Gclose(geom_group);
    free(geometry.height);
    return -1;
  }
  H5Dclose(dset_n_ij);

  // Write delta as scalar float64.
  hid_t dset_delta = H5Dcreate2(geom_group, "delta", H5T_IEEE_F64LE, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset_delta < 0 || H5Dwrite(dset_delta, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &geometry.delta) < 0) {
    fprintf(stderr, "Error writing geometry dataset delta\n");
    if (dset_delta >= 0)
      H5Dclose(dset_delta);
    H5Sclose(scalar_space);
    H5Gclose(geom_group);
    free(geometry.height);
    return -1;
  }
  H5Dclose(dset_delta);
  H5Sclose(scalar_space);

  // Write height as 1D float32 array.
  float *height_f32 = (float *)malloc(N_z * sizeof(float));
  if (height_f32 == NULL) {
    fprintf(stderr, "Error allocating memory for float32 height array\n");
    H5Gclose(geom_group);
    free(geometry.height);
    return -1;
  }

  for (int k = 0; k < N_z; ++k) {
    height_f32[k] = (float)geometry.height[k];
  }

  hsize_t height_dims[1] = {(hsize_t)N_z};
  hid_t   height_space   = H5Screate_simple(1, height_dims, NULL);
  if (height_space < 0) {
    fprintf(stderr, "Error creating dataspace for height dataset\n");
    free(height_f32);
    H5Gclose(geom_group);
    free(geometry.height);
    return -1;
  }

  hid_t dset_height =
      H5Dcreate2(geom_group, "height", H5T_IEEE_F32LE, height_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset_height < 0 || H5Dwrite(dset_height, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, height_f32) < 0) {
    fprintf(stderr, "Error writing geometry dataset height\n");
    if (dset_height >= 0)
      H5Dclose(dset_height);
    H5Sclose(height_space);
    free(height_f32);
    H5Gclose(geom_group);
    free(geometry.height);
    return -1;
  }
  H5Dclose(dset_height);
  H5Sclose(height_space);
  free(height_f32);
  // H5LTmake_dataset_double(geom_group, "delta_j", 1, (hsize_t[]){1}, &geometry.delta_j);

  H5Gclose(geom_group);
  free(geometry.height);

  return 0;
}

//============================================================================
// read_geometry_from_hdf5
//============================================================================
int
THDF_read_geometry_from_hdf5(hid_t    file,  //
                             int     *N_x,   //
                             int     *N_y,   //
                             int     *N_z,   //
                             double **heights) {
  hid_t geom_group = H5Gopen2(file, "/geometry", H5P_DEFAULT);
  if (geom_group < 0) {
    fprintf(stderr, "Error opening geometry group\n");
    return -1;
  }
  H5LTread_dataset_int(geom_group, "N_ij", N_x);
  *N_y = *N_x;
  H5LTread_dataset_int(geom_group, "N_height", N_z);
  *heights = (double *)malloc((*N_z) * sizeof(double));
  if (*heights == NULL) {
    fprintf(stderr, "Error allocating memory for height array\n");
    H5Gclose(geom_group);
    return EXIT_FAILURE;
  }
  H5LTread_dataset_double(geom_group, "height", *heights);
  H5Gclose(geom_group);
  return EXIT_SUCCESS;
}

//============================================================================
// write_frequency_grid_to_hdf5
//============================================================================
int
write_frequency_grid_to_hdf5(hid_t file, int N_frequencies, double *frequencies) {
  THDF_frequency_grid_t freq_grid = {N_frequencies, NULL, N_frequencies / 2};

  freq_grid.frequency = (double *)malloc(N_frequencies * sizeof(double));
  if (freq_grid.frequency == NULL) {
    fprintf(stderr, "Error allocating memory for frequency array\n");
    return EXIT_FAILURE;
  }

  for (int i = 0; i < N_frequencies; i++) {
    freq_grid.frequency[i] = frequencies[i];
  }

  // Write frequency grid data to HDF5 file
  hid_t freq_group = H5Gcreate2(file, "/frequency_grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (freq_group < 0) {
    fprintf(stderr, "Error creating frequency grid group\n");
    free(freq_grid.frequency);
    return EXIT_FAILURE;
  }

  H5LTmake_dataset_double(freq_group, "frequency", 1, (hsize_t[]){N_frequencies}, freq_grid.frequency);
  H5LTmake_dataset_int(freq_group, "N_frequencies", 1, (hsize_t[]){1}, &freq_grid.N_frequencies);

  H5Gclose(freq_group);
  free(freq_grid.frequency);

  return EXIT_SUCCESS;
}

//============================================================================
// read_frequency_grid_from_hdf5
//============================================================================
int
THDF_read_frequency_grid_from_hdf5(hid_t file, int *N_frequencies, double **frequencies) {

  hid_t freq_group = H5Gopen2(file, "/frequency_grid", H5P_DEFAULT);
  if (freq_group < 0) {
    fprintf(stderr, "Error opening frequency grid group\n");
    return EXIT_FAILURE;
  }
  H5LTread_dataset_int(freq_group, "N_frequencies", N_frequencies);

  *frequencies = (double *)malloc((*N_frequencies) * sizeof(double));
  if (*frequencies == NULL) {
    fprintf(stderr, "Error the memory for frequency array is not allocated\n");
    H5Gclose(freq_group);
    return EXIT_FAILURE;
  }
  H5LTread_dataset_double(freq_group, "frequency", *frequencies);
  H5Gclose(freq_group);
  return EXIT_SUCCESS;
}

//============================================================================
// create_atom_compound_type
//============================================================================
hid_t
create_atom_compound_type(void) {
  hid_t atom_type = H5Tcreate(H5T_COMPOUND, sizeof(THDF_atom_two_levels_t));
  H5Tinsert(atom_type, "atomic_number", HOFFSET(THDF_atom_two_levels_t, atomic_number), H5T_NATIVE_INT);
  H5Tinsert(atom_type, "atomic_mass", HOFFSET(THDF_atom_two_levels_t, atomic_mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(atom_type, "E_lower", HOFFSET(THDF_atom_two_levels_t, E_lower), H5T_NATIVE_DOUBLE);
  H5Tinsert(atom_type, "E_upper", HOFFSET(THDF_atom_two_levels_t, E_upper), H5T_NATIVE_DOUBLE);
  H5Tinsert(atom_type, "g_lower", HOFFSET(THDF_atom_two_levels_t, g_lower), H5T_NATIVE_DOUBLE);
  H5Tinsert(atom_type, "g_upper", HOFFSET(THDF_atom_two_levels_t, g_upper), H5T_NATIVE_DOUBLE);
  H5Tinsert(atom_type, "jl2", HOFFSET(THDF_atom_two_levels_t, jl2), H5T_NATIVE_INT);
  H5Tinsert(atom_type, "ju2", HOFFSET(THDF_atom_two_levels_t, ju2), H5T_NATIVE_INT);
  H5Tinsert(atom_type, "Aul", HOFFSET(THDF_atom_two_levels_t, Aul), H5T_NATIVE_DOUBLE);
  H5Tinsert(atom_type, "a_coef_D2", HOFFSET(THDF_atom_two_levels_t, a_coef_D2), H5T_NATIVE_DOUBLE);
  H5Tinsert(atom_type, "b_coef_D2", HOFFSET(THDF_atom_two_levels_t, b_coef_D2), H5T_NATIVE_DOUBLE);
  return atom_type;
}

//============================================================================
// write_atom_to_hdf5
//============================================================================
int
write_atom_to_hdf5(hid_t file, THDF_atom_two_levels_t *atom) {
  hsize_t dims[1]   = {1};
  hid_t   atom_type = create_atom_compound_type();
  hid_t   dataspace = H5Screate_simple(1, dims, NULL);
  hid_t   dataset   = H5Dcreate2(file, "/atom_data", atom_type, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0) {
    fprintf(stderr, "Error creating atom_data dataset\n");
    H5Sclose(dataspace);
    H5Tclose(atom_type);
    return EXIT_FAILURE;
  }
  H5Dwrite(dataset, atom_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, atom);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Tclose(atom_type);
  return EXIT_SUCCESS;
}

//============================================================================
// read_atom_from_hdf5
//============================================================================
int
THDF_read_atom_from_hdf5(hid_t file, THDF_atom_two_levels_t *atom) {
  hid_t dataset = H5Dopen2(file, "/atom_data", H5P_DEFAULT);
  if (dataset < 0) {
    fprintf(stderr, "Error opening atom_data dataset\n");
    return EXIT_FAILURE;
  }
  hid_t  atom_type = create_atom_compound_type();
  herr_t status    = H5Dread(dataset, atom_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, atom);
  H5Tclose(atom_type);
  H5Dclose(dataset);
  if (status < 0) {
    fprintf(stderr, "Error reading atom_data dataset\n");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

//============================================================================
// create_atom_two_terms_compound_type
//============================================================================
hid_t
create_atom_two_terms_compound_type(void) {
  hsize_t E_dims[1]    = {MAX_ATOM_TERMS_PLACEHOLDER};
  hid_t   E_array_type = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, E_dims);

  hid_t atom_type = H5Tcreate(H5T_COMPOUND, sizeof(THDF_atom_two_terms_t));
  H5Tinsert(atom_type, "atomic_number", HOFFSET(THDF_atom_two_terms_t, atomic_number), H5T_NATIVE_INT);
  H5Tinsert(atom_type, "atomic_mass", HOFFSET(THDF_atom_two_terms_t, atomic_mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(atom_type, "Aul", HOFFSET(THDF_atom_two_terms_t, Aul), H5T_NATIVE_DOUBLE);
  H5Tinsert(atom_type, "L2_upper", HOFFSET(THDF_atom_two_terms_t, L2_upper), H5T_NATIVE_INT);
  H5Tinsert(atom_type, "L2_lower", HOFFSET(THDF_atom_two_terms_t, L2_lower), H5T_NATIVE_INT);
  H5Tinsert(atom_type, "S2", HOFFSET(THDF_atom_two_terms_t, S2), H5T_NATIVE_INT);
  H5Tinsert(atom_type, "g_lower", HOFFSET(THDF_atom_two_terms_t, g_lower), H5T_NATIVE_DOUBLE);
  H5Tinsert(atom_type, "g_upper", HOFFSET(THDF_atom_two_terms_t, g_upper), H5T_NATIVE_DOUBLE);
  H5Tinsert(atom_type, "E_size", HOFFSET(THDF_atom_two_terms_t, E_size), H5T_NATIVE_INT);
  H5Tinsert(atom_type, "E_lower", HOFFSET(THDF_atom_two_terms_t, E_lower), E_array_type);
  H5Tinsert(atom_type, "E_upper", HOFFSET(THDF_atom_two_terms_t, E_upper), E_array_type);

  H5Tclose(E_array_type);
  return atom_type;
}

//============================================================================
// write_atom_two_terms_to_hdf5
//============================================================================
int
write_atom_two_terms_to_hdf5(hid_t file, THDF_atom_two_terms_t *atom) {
  if (atom->E_size < 0 || atom->E_size > MAX_ATOM_TERMS_PLACEHOLDER) {
    fprintf(stderr, "Error: E_size (%d) out of range [0, %d]\n", atom->E_size, MAX_ATOM_TERMS_PLACEHOLDER);
    return EXIT_FAILURE;
  }

  // Levels beyond E_size are unused: zero them out before writing so no
  // uninitialized/stale values are persisted to the file.
  for (int i = atom->E_size; i < MAX_ATOM_TERMS_PLACEHOLDER; i++) {
    atom->E_lower[i] = 0.0;
    atom->E_upper[i] = 0.0;
  }

  hsize_t dims[1]   = {1};
  hid_t   atom_type = create_atom_two_terms_compound_type();
  hid_t   dataspace = H5Screate_simple(1, dims, NULL);
  hid_t   dataset = H5Dcreate2(file, "/atom_two_terms_data", atom_type, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0) {
    fprintf(stderr, "Error creating atom_two_terms_data dataset\n");
    H5Sclose(dataspace);
    H5Tclose(atom_type);
    return EXIT_FAILURE;
  }
  H5Dwrite(dataset, atom_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, atom);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Tclose(atom_type);
  return EXIT_SUCCESS;
}

//============================================================================
// THDF_read_atom_two_terms_from_hdf5
//============================================================================
int
THDF_read_atom_two_terms_from_hdf5(hid_t file, THDF_atom_two_terms_t *atom) {
  hid_t dataset = H5Dopen2(file, "/atom_two_terms_data", H5P_DEFAULT);
  if (dataset < 0) {
    fprintf(stderr, "Error opening atom_two_terms_data dataset\n");
    return EXIT_FAILURE;
  }
  hid_t  atom_type = create_atom_two_terms_compound_type();
  herr_t status    = H5Dread(dataset, atom_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, atom);
  H5Tclose(atom_type);
  H5Dclose(dataset);
  if (status < 0) {
    fprintf(stderr, "Error reading atom_two_terms_data dataset\n");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

//============================================================================
// THDF_is_two_term_atom
//============================================================================
bool
THDF_is_two_term_atom(hid_t file) {
  htri_t exists = H5Lexists(file, "/atom_two_terms_data", H5P_DEFAULT);
  return exists > 0;
}

//============================================================================
// write_atmos_3D_data_to_hdf5
//============================================================================
int
write_atmos_3D_data_to_hdf5(hid_t file, THDF_atmos_t *atmos_data, int N_x, int N_y, int N_z) {
  hsize_t dims[3] = {N_x, N_y, N_z};

  // Create the compound datatype for hdf_atmos
  hid_t atmos_type = create_hdf_atmos_compound_type();
  // Create the dataset at root
  hid_t dataspace = H5Screate_simple(3, dims, NULL);
  hid_t dataset   = H5Dcreate2(file, "atmos", atmos_type, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0) {
    fprintf(stderr, "Error creating atmosphere dataset\n");
    H5Tclose(atmos_type);
    return EXIT_FAILURE;
  }
  H5Dwrite(dataset, atmos_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, atmos_data);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Tclose(atmos_type);
  return EXIT_SUCCESS;
}

//============================================================================
// read_atmos_3D_data_from_hdf5
//============================================================================
int
THDF_read_atmos_3D_data_from_hdf5(hid_t file, THDF_atmos_t **atmos_data, int *N_x, int *N_y, int *N_z) {
  // Open the dataset at root
  hid_t dataset = H5Dopen2(file, "atmos", H5P_DEFAULT);
  if (dataset < 0) {
    fprintf(stderr, "Error opening atmosphere dataset\n");
    return EXIT_FAILURE;
  }

  // Get dataspace and read dimensions
  hid_t   dataspace = H5Dget_space(dataset);
  hsize_t dims[3];
  H5Sget_simple_extent_dims(dataspace, dims, NULL);
  H5Sclose(dataspace);

  *N_x = (int)dims[0];
  *N_y = (int)dims[1];
  *N_z = (int)dims[2];

  size_t total_size = (*N_x) * (*N_y) * (*N_z);
  *atmos_data       = (THDF_atmos_t *)malloc(total_size * sizeof(THDF_atmos_t));
  if (*atmos_data == NULL) {
    fprintf(stderr, "Error allocating memory for atmosphere data\n");
    H5Dclose(dataset);
    return EXIT_FAILURE;
  }

  // Create the compound datatype for hdf_atmos
  hid_t atmos_type = create_hdf_atmos_compound_type();

  // Read the data
  H5Dread(dataset, atmos_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, *atmos_data);
  H5Dclose(dataset);
  H5Tclose(atmos_type);
  return EXIT_SUCCESS;
}

//============================================================================
// read_atmos_3D_data_from_hdf5 (partial read)
//============================================================================
int
THDF_read_atmos_3D_block_hdf5(hid_t file, THDF_atmos_t **atmos_data,    //
                              int start_i, int start_j, int start_k,    //
                              int count_i, int count_j, int count_k) {  //

  // Open the dataset at root
  hid_t dataset = H5Dopen2(file, "atmos", H5P_DEFAULT);
  if (dataset < 0) {
    fprintf(stderr, "Error opening atmosphere dataset\n");
    return EXIT_FAILURE;
  }

  hsize_t start[3] = {start_i, start_j, start_k};
  hsize_t count[3] = {count_i, count_j, count_k};

  size_t total_size = count_i * count_j * count_k;
  *atmos_data       = (THDF_atmos_t *)malloc(total_size * sizeof(THDF_atmos_t));
  if (*atmos_data == NULL) {
    fprintf(stderr, "Error allocating memory for atmosphere data\n");
    H5Dclose(dataset);
    return EXIT_FAILURE;
  }

  // Create the compound datatype for hdf_atmos
  hid_t atmos_type = create_hdf_atmos_compound_type();

  // Create memory space for the data to be read
  hid_t memspace = H5Screate_simple(3, count, NULL);

  // Select hyperslab in file dataspace
  hid_t filespace = H5Dget_space(dataset);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);

  // Read the data
  H5Dread(dataset, atmos_type, memspace, filespace, H5P_DEFAULT, *atmos_data);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dataset);
  H5Tclose(atmos_type);
  return EXIT_SUCCESS;
}

//============================================================================
// read_atmos_3D_data_from_hdf5 (partial read)
//============================================================================
int
THDF_read_atmos_3D_point_hdf5(hid_t file, THDF_atmos_t *atmos_data,     //
                              int start_i, int start_j, int start_k) {  //

  // Open the dataset at root
  hid_t dataset = H5Dopen2(file, "atmos", H5P_DEFAULT);
  if (dataset < 0) {
    fprintf(stderr, "Error opening atmosphere dataset\n");
    return EXIT_FAILURE;
  }

  hsize_t start[3] = {start_i, start_j, start_k};
  hsize_t count[3] = {1, 1, 1};

  // Create the compound datatype for hdf_atmos
  hid_t atmos_type = create_hdf_atmos_compound_type();

  // Create memory space for the data to be read
  hid_t memspace = H5Screate_simple(3, count, NULL);

  // Select hyperslab in file dataspace
  hid_t filespace = H5Dget_space(dataset);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);

  // Read the data
  H5Dread(dataset, atmos_type, memspace, filespace, H5P_DEFAULT, atmos_data);

  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Tclose(atmos_type);
  H5Dclose(dataset);
  return EXIT_SUCCESS;
}

//============================================================================
// create_atmosphere_dataset
//============================================================================
THDF_atmos_dataset_handler_t *
THDF_create_atmosphere_dataset(hid_t file, int N_x, int N_y, int N_z) {
  THDF_atmos_dataset_handler_t *atmos_dset = (THDF_atmos_dataset_handler_t *)malloc(sizeof(THDF_atmos_dataset_handler_t));
  if (!atmos_dset) {
    fprintf(stderr, "Error allocating memory for atmosphere dataset struct\n");
    return NULL;
  }

  // Initialize struct
  atmos_dset->file_id = file;
  atmos_dset->N_x     = N_x;
  atmos_dset->N_y     = N_y;
  atmos_dset->N_z     = N_z;
  atmos_dset->is_open = 0;

  hsize_t dims[3] = {N_x, N_y, N_z};

  atmos_dset->group_id = -1;  // dataset lives at root, no group

  // Create the compound datatype for hdf_atmos
  atmos_dset->datatype_id = create_hdf_atmos_compound_type();

  // Create dataspace and dataset at root
  atmos_dset->dataspace_id = H5Screate_simple(3, dims, NULL);
  atmos_dset->dataset_id =
      H5Dcreate2(file, "atmos", atmos_dset->datatype_id, atmos_dset->dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (atmos_dset->dataset_id < 0) {
    fprintf(stderr, "Error creating atmosphere dataset\n");
    H5Tclose(atmos_dset->datatype_id);
    H5Sclose(atmos_dset->dataspace_id);
    free(atmos_dset);
    return NULL;
  }

  atmos_dset->is_open = 1;
  return atmos_dset;
}

//============================================================================
// write_atmosphere_data_to_dataset
//============================================================================
int
THDF_write_atmosphere_data_to_dataset(THDF_atmos_dataset_handler_t *atmos_dset, THDF_atmos_t *atmos_data, hsize_t start_i,
                                      hsize_t start_j, hsize_t start_k, hsize_t count_i, hsize_t count_j,
                                      hsize_t count_k) {
  if (!atmos_dset || !atmos_dset->is_open) {
    fprintf(stderr, "Error: atmosphere dataset is not open\n");
    return -1;
  }

  // Define hyperslab for partial writing
  hsize_t start[3] = {start_i, start_j, start_k};
  hsize_t count[3] = {count_i, count_j, count_k};

  hid_t memspace = H5Screate_simple(3, count, NULL);
  H5Sselect_hyperslab(atmos_dset->dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

  herr_t status = H5Dwrite(atmos_dset->dataset_id, atmos_dset->datatype_id, memspace, atmos_dset->dataspace_id,
                           H5P_DEFAULT, atmos_data);

  H5Sclose(memspace);

  if (status < 0) {
    fprintf(stderr, "Error writing atmosphere data to dataset\n");
    return -1;
  }

  return 0;
}

//============================================================================
// close_atmosphere_dataset
//============================================================================
int
THDF_close_atmosphere_dataset(THDF_atmos_dataset_handler_t *atmos_dset) {
  if (!atmos_dset) {
    return -1;
  }

  if (atmos_dset->is_open) {
    H5Dclose(atmos_dset->dataset_id);
    H5Sclose(atmos_dset->dataspace_id);
    H5Tclose(atmos_dset->datatype_id);
    if (atmos_dset->group_id >= 0)
      H5Gclose(atmos_dset->group_id);
    atmos_dset->is_open = 0;
  }

  free(atmos_dset);
  return 0;
}

//============================================================================
// read_main_demo_from_file
//============================================================================
// read_main_demo_from_file
//============================================================================
int
THDF_read_main_demo_from_file(const char *filename, const bool verbose) {

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

  int     N_frequencies;
  double *frequencies = NULL;
  if (THDF_read_frequency_grid_from_hdf5(file, &N_frequencies, &frequencies) < 0) {
    fprintf(stderr, "Error reading frequency grid data\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }

  printf("Frequency grid: N_frequencies=%d\n", N_frequencies);
  printf("Frequencies:\n");
  if (verbose) {
    for (int i = 0; i < N_frequencies; i++) {
      printf("  frequency[%d] = %e Hz\n", i, frequencies[i]);
    }
  }
  free(frequencies);

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

  // THDF_read_atmos_3D_data_from_hdf5(file, &atmos_data, &N_x, &N_y, &N_z);

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
  THDF_print_atmos_point(&atmos_data[0]);

  free(atmos_data);
  H5Fclose(file);

  return EXIT_SUCCESS;
}