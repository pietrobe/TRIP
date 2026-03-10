#include "hdf_output_field.h"

#include <stdio.h>
#include <stdlib.h>

//============================================================================
// MPI Utility Functions
//============================================================================

hid_t
THDF_open_file_MPI(const char *filename, MPI_Comm mpi_comm) {
  hid_t file_id, plist_id;

  // Create file access property list
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  if (plist_id < 0) {
    fprintf(stderr, "Error creating file access property list\n");
    return -1;
  }

  // Set MPI-IO file access
  herr_t status = H5Pset_fapl_mpio(plist_id, mpi_comm, MPI_INFO_NULL);
  if (status < 0) {
    fprintf(stderr, "Error setting MPI-IO file access property\n");
    H5Pclose(plist_id);
    return -1;
  }

  // Open the file
  file_id = H5Fopen(filename, H5F_ACC_RDWR, plist_id);
  if (file_id < 0) {
    fprintf(stderr, "Error opening file '%s' in MPI mode\n", filename);
    H5Pclose(plist_id);
    return -1;
  }

  // Close the property list
  H5Pclose(plist_id);

  return file_id;
}

//============================================================================
// File Creation Functions
//============================================================================

hid_t
THDF_open_file(const char *filename) {  //
  hid_t file_id, plist_id;

  // Create file access property list
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  if (plist_id < 0) {
    fprintf(stderr, "Error creating file access property list\n");
    return -1;
  }

  // Create the file (truncate if exists)
  file_id = H5Fcreate(filename,       //
                      H5F_ACC_TRUNC,  //
                      H5P_DEFAULT,    //
                      plist_id);      //
  if (file_id < 0) {
    fprintf(stderr, "Error creating file '%s'\n", filename);
    H5Pclose(plist_id);
    return -1;
  }

  // Close the property list
  H5Pclose(plist_id);

  return file_id;
}

void                              //
THDF_close_file(hid_t file_id) {  //
  H5Fclose(file_id);
}

//============================================================================
// Type and Datatype Functions
//============================================================================

THDF_float_type_t
THDF_get_float_type(void) {
  switch (sizeof(THDF_float_t)) {
    case 8:
      return HDF_OUT_FLOAT64;
      break;

    case 4:
      return HDF_OUT_FLOAT32;
      break;

    default:
      return HDF_OUT_NONE;
      break;
  }
}

hid_t
THDF_get_hdf_float_datatype(void) {
  switch (sizeof(THDF_float_t)) {
    case 8:
      return H5T_NATIVE_DOUBLE;
      break;

    case 4:
      return H5T_NATIVE_FLOAT;
      break;

    default:
      return -1;
      break;
  }
}

hid_t
THDF_get_hdf_float32_datatype(void) {
  return H5T_NATIVE_FLOAT;
}

THDF_n_float_type_t
THDF_get_n_float_type(void) {
  switch (sizeof(THDF_n_float_t)) {
    case 8:
      return HDF_OUT_NORM_FLOAT64;
      break;

    case 4:
      return HDF_OUT_NORM_FLOAT32;
      break;

    default:
      return HDF_OUT_NORM_NONE;
      break;
  }
}

hid_t
THDF_get_hdf_n_float_datatype(void) {
  switch (sizeof(THDF_n_float_t)) {
    case 8:
      return H5T_NATIVE_DOUBLE;
      break;

    case 4:
      return H5T_NATIVE_FLOAT;
      break;

    default:
      return -1;
      break;
  }
}

//============================================================================
// THDF_create_field_handler_mpi
//============================================================================
/////////////////////////////////////////////////////////////
/// THDF_create_field_handler_mpi
/////////////////////////////////////////////////////////////
THDF_field_handler_t *
THDF_create_field_handler_mpi(hid_t file,           //
                              int   N_x,            //
                              int   N_y,            //
                              int   N_incl,         //
                              int   N_azimuth,      //
                              int   N_frequencies) {  //
                                                    //
  THDF_field_handler_t *output_dset = (THDF_field_handler_t *)malloc(sizeof(THDF_field_handler_t));

  if (!output_dset) {
    fprintf(stderr, "Error allocating memory for output field dataset struct\n");
    return NULL;
  }

  // Initialize struct
  output_dset->file_id       = file;
  output_dset->N_x           = N_x;
  output_dset->N_y           = N_y;
  output_dset->N_incl        = N_incl;
  output_dset->N_azimuth     = N_azimuth;
  output_dset->N_frequencies = N_frequencies;
  output_dset->is_open       = 0;

  // const int N_stokes = 4;  // I, Q, U, V

  hsize_t dims[5] = {N_x, N_y, N_incl, N_azimuth, N_frequencies};

  // Check if output_field group already exists, if so, open it
  H5E_BEGIN_TRY { output_dset->group_id = H5Gopen2(file, "/emergent_field", H5P_DEFAULT); }
  H5E_END_TRY;

  // If group doesn't exist, create it (collective operation)
  if (output_dset->group_id < 0) {
    output_dset->group_id = H5Gcreate2(file, "/emergent_field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (output_dset->group_id < 0) {
      fprintf(stderr, "Error creating output_field group with MPI\n");
      free(output_dset);
      return NULL;
    }
  }

  // Create the datatype for double
  output_dset->datatype_id = H5Tcopy(THDF_get_hdf_float_datatype());
  // Create dataspace and dataset with MPI support (all processes collectively)
  output_dset->dataspace_id_I = H5Screate_simple(5, dims, NULL);
  output_dset->dataset_id_I   = H5Dcreate2(output_dset->group_id, "emergent_I", output_dset->datatype_id,
                                           output_dset->dataspace_id_I, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  output_dset->dataspace_id_Q = H5Screate_simple(5, dims, NULL);
  output_dset->dataset_id_Q   = H5Dcreate2(output_dset->group_id, "emergent_QI_pc", output_dset->datatype_id,
                                           output_dset->dataspace_id_Q, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  output_dset->dataspace_id_U = H5Screate_simple(5, dims, NULL);
  output_dset->dataset_id_U   = H5Dcreate2(output_dset->group_id, "emergent_UI_pc", output_dset->datatype_id,
                                           output_dset->dataspace_id_U, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  output_dset->dataspace_id_V = H5Screate_simple(5, dims, NULL);
  output_dset->dataset_id_V   = H5Dcreate2(output_dset->group_id, "emergent_VI_pc", output_dset->datatype_id,
                                           output_dset->dataspace_id_V, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (output_dset->dataset_id_I < 0 || output_dset->dataset_id_Q < 0 || output_dset->dataset_id_U < 0 ||
      output_dset->dataset_id_V < 0) {
    fprintf(stderr, "Error creating output field dataset with MPI\n");
    H5Tclose(output_dset->datatype_id);
    H5Sclose(output_dset->dataspace_id_I);
    H5Gclose(output_dset->group_id);
    free(output_dset);
    return NULL;
  }
  output_dset->is_open = 1;
  return output_dset;
}

/////////////////////////////////////////////////////////////
/// THDF_close_field_handler_mpi
/////////////////////////////////////////////////////////////
void
THDF_close_field_handler_mpi(THDF_field_handler_t *output_dset) {
  if (!output_dset) {
    return;
  }

  if (output_dset->is_open) {
    H5Dclose(output_dset->dataset_id_I);
    H5Dclose(output_dset->dataset_id_Q);
    H5Dclose(output_dset->dataset_id_U);
    H5Dclose(output_dset->dataset_id_V);
    H5Sclose(output_dset->dataspace_id_I);
    H5Sclose(output_dset->dataspace_id_Q);
    H5Sclose(output_dset->dataspace_id_U);
    H5Sclose(output_dset->dataspace_id_V);

    H5Tclose(output_dset->datatype_id);
    H5Gclose(output_dset->group_id);
    output_dset->is_open = 0;
  }

  free(output_dset);
}

/////////////////////////////////////////////////////////////
/// THDF_write_field_dataset_to_hdf5
/////////////////////////////////////////////////////////////
int
THDF_write_field_dataset_to_hdf5(THDF_field_handler_t *output_dset,    //
                                 THDF_field_t         *output_field,   //
                                 hsize_t               start_i,        //
                                 hsize_t               start_j,        //
                                 hsize_t               start_incl,     //
                                 hsize_t               start_azimuth,  //
                                 hsize_t               count_i,        //
                                 hsize_t               count_j,        //
                                 hsize_t               count_incl,     //
                                 hsize_t               count_azimuth,  //
                                 hsize_t               count_frequencies) {
  if (!output_dset || !output_dset->is_open) {
    fprintf(stderr, "Error: output field dataset is not open\n");
    return -1;
  }

  // 5D hyperslab: [i, j, stokes, angular_index, frequencies]
  // Each dataset is separate, so stokes index is always 0
  // hsize_t angular_index = start_incl * output_dset->N_azimuth + start_azimuth;

  hsize_t start[5] = {start_i, start_j, start_incl, start_azimuth, 0};
  hsize_t count[5] = {count_i, count_j, count_incl, count_azimuth, count_frequencies};

  // Memory space should match data layout: [x, y, incl, azimuth, frequencies]
  hid_t memspace = H5Screate_simple(5, count, NULL);

  // Create separate dataspaces for each Stokes parameter
  hid_t dataspace_I = H5Dget_space(output_dset->dataset_id_I);
  hid_t dataspace_Q = H5Dget_space(output_dset->dataset_id_Q);
  hid_t dataspace_U = H5Dget_space(output_dset->dataset_id_U);
  hid_t dataspace_V = H5Dget_space(output_dset->dataset_id_V);

  H5Sselect_hyperslab(dataspace_I, H5S_SELECT_SET, start, NULL, count, NULL);
  H5Sselect_hyperslab(dataspace_Q, H5S_SELECT_SET, start, NULL, count, NULL);
  H5Sselect_hyperslab(dataspace_U, H5S_SELECT_SET, start, NULL, count, NULL);
  H5Sselect_hyperslab(dataspace_V, H5S_SELECT_SET, start, NULL, count, NULL);

  // Write each Stokes parameter separately
  herr_t status_I =                        //
      H5Dwrite(output_dset->dataset_id_I,  //
               output_dset->datatype_id,   //
               memspace,                   //
               dataspace_I,                //
               H5P_DEFAULT,                //
               output_field->stokes_I);    //
  HDF5_CHECK_WRITE(status_I, "Stokes I", memspace);

  herr_t status_Q =                        //
      H5Dwrite(output_dset->dataset_id_Q,  //
               output_dset->datatype_id,   //
               memspace,                   //
               dataspace_Q,                //
               H5P_DEFAULT,                //
               output_field->stokes_QI);   //
  HDF5_CHECK_WRITE(status_Q, "Stokes QI", memspace);
  herr_t status_U =                        //
      H5Dwrite(output_dset->dataset_id_U,  //
               output_dset->datatype_id,   //
               memspace,                   //
               dataspace_U,                //
               H5P_DEFAULT,                //
               output_field->stokes_UI);   //
  HDF5_CHECK_WRITE(status_U, "Stokes UI", memspace);

  herr_t status_V =                        //
      H5Dwrite(output_dset->dataset_id_V,  //
               output_dset->datatype_id,   //
               memspace,                   //
               dataspace_V,                //
               H5P_DEFAULT,                //
               output_field->stokes_VI);   //
  HDF5_CHECK_WRITE(status_V, "Stokes VI", memspace);

  H5Sclose(memspace);
  H5Sclose(dataspace_I);
  H5Sclose(dataspace_Q);
  H5Sclose(dataspace_U);
  H5Sclose(dataspace_V);
  return 0;
}

/////////////////////////////////////////////////////////////
/// THDF_read_field_from_hdf5
/////////////////////////////////////////////////////////////
THDF_field_t
THDF_read_field_from_hdf5(hid_t     file,           //
                          const int N_frequencies,  //
                          const int x_i,            //
                          const int y_i,            //
                          const int incl_i,         //
                          const int azimuth_i) {    //

  THDF_field_t output_field = THDF_make_empty_field();

  hid_t group_id = H5Gopen2(file, "/emergent_field", H5P_DEFAULT);
  if (group_id < 0) {
    fprintf(stderr, "Error opening output_field group\n");
    return output_field;
  }

  hid_t dataset_id_I = H5Dopen2(group_id, "emergent_I", H5P_DEFAULT);
  hid_t dataset_id_Q = H5Dopen2(group_id, "emergent_QI_pc", H5P_DEFAULT);
  hid_t dataset_id_U = H5Dopen2(group_id, "emergent_UI_pc", H5P_DEFAULT);
  hid_t dataset_id_V = H5Dopen2(group_id, "emergent_VI_pc", H5P_DEFAULT);

  if (dataset_id_I < 0 || dataset_id_Q < 0 || dataset_id_U < 0 || dataset_id_V < 0) {
    fprintf(stderr, "Error opening output field datasets\n");
    fprintf(stderr, "dataset_id_I: %ld, dataset_id_Q: %ld, dataset_id_U: %ld, dataset_id_V: %ld\n", dataset_id_I,
            dataset_id_Q, dataset_id_U, dataset_id_V);
    H5Gclose(group_id);
    return output_field;
  }

  // Get dataspace to define hyperslab
  hid_t dataspace_id_I = H5Dget_space(dataset_id_I);
  if (dataspace_id_I < 0) {
    fprintf(stderr, "Error getting dataspace for Stokes I\n");
    H5Dclose(dataset_id_I);
    H5Dclose(dataset_id_Q);
    H5Dclose(dataset_id_U);
    H5Dclose(dataset_id_V);
    H5Gclose(group_id);
    return output_field;
  }

  output_field.index_i       = x_i;
  output_field.index_j       = y_i;
  output_field.index_incl    = incl_i;
  output_field.index_azimuth = azimuth_i;
  output_field.stokes_I      = (THDF_float_t *)malloc(N_frequencies * sizeof(THDF_float_t));
  output_field.stokes_QI     = (THDF_float_t *)malloc(N_frequencies * sizeof(THDF_float_t));
  output_field.stokes_UI     = (THDF_float_t *)malloc(N_frequencies * sizeof(THDF_float_t));
  output_field.stokes_VI     = (THDF_float_t *)malloc(N_frequencies * sizeof(THDF_float_t));

  hid_t HDF_DATATYPE_ID = THDF_get_hdf_float_datatype();

  // Define hyperslab to read
  hsize_t start[5] = {x_i, y_i, incl_i, azimuth_i, 0};
  hsize_t count[5] = {1, 1, 1, 1, N_frequencies};
  H5Sselect_hyperslab(dataspace_id_I, H5S_SELECT_SET, start, NULL, count, NULL);
  hid_t memspace = H5Screate_simple(2, (hsize_t[]){1, N_frequencies}, NULL);
  // Read Stokes I
  herr_t status = H5Dread(dataset_id_I, HDF_DATATYPE_ID, memspace, dataspace_id_I, H5P_DEFAULT, output_field.stokes_I);
  HDF5_CHECK_WRITE_OF(status, "Stokes I", memspace, output_field);

  // Read Stokes Q
  hid_t dataspace_id_Q = H5Dget_space(dataset_id_Q);
  H5Sselect_hyperslab(dataspace_id_Q, H5S_SELECT_SET, start, NULL, count, NULL);
  status = H5Dread(dataset_id_Q, HDF_DATATYPE_ID, memspace, dataspace_id_Q, H5P_DEFAULT, output_field.stokes_QI);
  HDF5_CHECK_WRITE_OF(status, "Stokes QI", memspace, output_field);

  // Read Stokes U
  hid_t dataspace_id_U = H5Dget_space(dataset_id_U);
  H5Sselect_hyperslab(dataspace_id_U, H5S_SELECT_SET, start, NULL, count, NULL);
  status = H5Dread(dataset_id_U, HDF_DATATYPE_ID, memspace, dataspace_id_U, H5P_DEFAULT, output_field.stokes_UI);
  HDF5_CHECK_WRITE_OF(status, "Stokes UI", memspace, output_field);

  // Read Stokes V
  hid_t dataspace_id_V = H5Dget_space(dataset_id_V);
  H5Sselect_hyperslab(dataspace_id_V, H5S_SELECT_SET, start, NULL, count, NULL);
  status = H5Dread(dataset_id_V, HDF_DATATYPE_ID, memspace, dataspace_id_V, H5P_DEFAULT, output_field.stokes_VI);
  HDF5_CHECK_WRITE_OF(status, "Stokes VI", memspace, output_field);

  H5Sclose(memspace);

  H5Sclose(dataspace_id_I);
  H5Sclose(dataspace_id_Q);
  H5Sclose(dataspace_id_U);
  H5Sclose(dataspace_id_V);

  H5Dclose(dataset_id_I);
  H5Dclose(dataset_id_Q);
  H5Dclose(dataset_id_U);
  H5Dclose(dataset_id_V);
  H5Gclose(group_id);
  return output_field;
}

/////////////////////////////////////////////////////////////
/// THDF_write_angular_grid_to_hdf5
/////////////////////////////////////////////////////////////
int                                                               //
THDF_write_angular_grid_to_hdf5(hid_t                      file,  //
                                const THDF_angular_grid_t *angular_grid) {

  hid_t angular_group = H5Gcreate2(file, "/emergent_angular_grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (angular_group < 0) {
    fprintf(stderr, "Error creating emergent angular grid group\n");
    return -1;
  }

  H5LTmake_dataset_int(angular_group, "N_directions", 1, (hsize_t[]){1},  //
                       &angular_grid->N_directions);

  H5LTmake_dataset_int(angular_group, "N_inclination_angles", 1, (hsize_t[]){1},  //
                       &angular_grid->N_inclination_angles);

  H5LTmake_dataset_int(angular_group, "N_azimuthal_angles", 1, (hsize_t[]){1},  //
                       &angular_grid->N_azimuthal_angles);

  H5LTmake_dataset_double(angular_group, "inclination_angles", 1, (hsize_t[]){angular_grid->N_inclination_angles},
                          angular_grid->inclination_angles);

  H5LTmake_dataset_double(angular_group, "azimuthal_angles", 1, (hsize_t[]){angular_grid->N_azimuthal_angles},
                          angular_grid->azimuthal_angles);

  H5LTmake_dataset_int(angular_group, "inclinations_indices", 1, (hsize_t[]){angular_grid->N_directions},
                       angular_grid->inclinations_indices);

  H5LTmake_dataset_int(angular_group, "azimuthal_indices", 1, (hsize_t[]){angular_grid->N_directions},
                       angular_grid->azimuthal_indices);
  H5Gclose(angular_group);
  return 0;
}  //

/////////////////////////////////////////////////////////////
/// THDF_read_angular_grid_from_hdf5
/////////////////////////////////////////////////////////////
THDF_angular_grid_t
THDF_read_angular_grid_from_hdf5(hid_t file) {
  THDF_angular_grid_t angular_grid = THDF_make_empty_angular_grid();

  hid_t angular_group = H5Gopen2(file, "/emergent_angular_grid", H5P_DEFAULT);
  if (angular_group < 0) {
    fprintf(stderr, "Error opening emergent angular grid group\n");
    return angular_grid;
  }

  H5LTread_dataset_int(angular_group, "N_directions", &angular_grid.N_directions);
  H5LTread_dataset_int(angular_group, "N_inclination_angles", &angular_grid.N_inclination_angles);
  H5LTread_dataset_int(angular_group, "N_azimuthal_angles", &angular_grid.N_azimuthal_angles);

  angular_grid.inclination_angles = (double *)malloc((angular_grid.N_directions) * sizeof(double));
  if (angular_grid.inclination_angles == NULL) {
    fprintf(stderr, "Error allocating memory for inclination angles array\n");
    H5Gclose(angular_group);
    return angular_grid;
  }
  H5LTread_dataset_double(angular_group, "inclination_angles", angular_grid.inclination_angles);

  angular_grid.azimuthal_angles = (double *)malloc((angular_grid.N_directions) * sizeof(double));
  if (angular_grid.azimuthal_angles == NULL) {
    fprintf(stderr, "Error allocating memory for azimuthal angles array\n");
    free(angular_grid.inclination_angles);
    H5Gclose(angular_group);
    return angular_grid;
  }
  H5LTread_dataset_double(angular_group, "azimuthal_angles", angular_grid.azimuthal_angles);

  angular_grid.inclinations_indices = (int *)malloc((angular_grid.N_directions) * sizeof(int));
  if (angular_grid.inclinations_indices == NULL) {
    fprintf(stderr, "Error allocating memory for inclinations indices array\n");
    free(angular_grid.inclination_angles);
    free(angular_grid.azimuthal_angles);
    H5Gclose(angular_group);
    return angular_grid;
  }
  H5LTread_dataset_int(angular_group, "inclinations_indices", angular_grid.inclinations_indices);

  angular_grid.azimuthal_indices = (int *)malloc((angular_grid.N_directions) * sizeof(int));
  if (angular_grid.azimuthal_indices == NULL) {
    fprintf(stderr, "Error allocating memory for azimuthal indices array\n");
    free(angular_grid.inclination_angles);
    free(angular_grid.azimuthal_angles);
    free(angular_grid.inclinations_indices);
    H5Gclose(angular_group);
    return angular_grid;
  }

  H5LTread_dataset_int(angular_group, "azimuthal_indices", angular_grid.azimuthal_indices);
  H5Gclose(angular_group);
  return angular_grid;
}

/////////////////////////////////////////////////////////////
/// THDF_write_frequencies_grid_to_hdf5
/////////////////////////////////////////////////////////////
int
THDF_write_frequencies_grid_to_hdf5(hid_t                          file,                //
                                    const THDF_frequencies_grid_t *frequencies_grid) {  //

  hid_t freq_group = H5Gcreate2(file, "/frequencies_grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (freq_group < 0) {
    fprintf(stderr, "Error creating output frequencies grid group\n");
    return -1;
  }
  H5LTmake_dataset_double(freq_group, "frequencies", 1, (hsize_t[]){frequencies_grid->N_frequencies},
                          frequencies_grid->frequencies);
  H5Gclose(freq_group);
  return 0;
}

///////////////////////////////////////////////////////////////
/// THDF_read_frequencies_grid_from_hdf5
/////////////////////////////////////////////////////////////
int
THDF_write_geometry_3D_to_hdf5(hid_t file, const THDF_geometry_3D_t *geometry) {
  hid_t geom_group = H5Gcreate2(file, "/geometry_3D", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (geom_group < 0) {
    fprintf(stderr, "Error creating geometry_3D group\n");
    return -1;
  }

  H5LTmake_dataset_int(geom_group, "N_x", 1, (hsize_t[]){1}, &geometry->N_x);
  H5LTmake_dataset_int(geom_group, "N_y", 1, (hsize_t[]){1}, &geometry->N_y);
  H5LTmake_dataset_int(geom_group, "N_z", 1, (hsize_t[]){1}, &geometry->N_z);
  H5LTmake_dataset_double(geom_group, "heights", 1, (hsize_t[]){(hsize_t)(geometry->N_z)}, geometry->heights);
  H5LTmake_dataset_double(geom_group, "delta", 1, (hsize_t[]){1}, &geometry->delta);

  H5Gclose(geom_group);
  return 0;
}

/////////////////////////////////////////////////////////////
/// THDF_read_frequencies_grid_from_hdf5
/////////////////////////////////////////////////////////////
THDF_frequencies_grid_t
THDF_read_frequencies_grid_from_hdf5(hid_t file) {
  THDF_frequencies_grid_t frequencies_grid = THDF_make_empty_frequencies_grid();

  hid_t freq_group = H5Gopen2(file, "/frequencies_grid", H5P_DEFAULT);
  if (freq_group < 0) {
    fprintf(stderr, "Error opening output frequencies grid group\n");
    return frequencies_grid;
  }
  H5LTread_dataset_int(freq_group, "N_frequencies", &frequencies_grid.N_frequencies);

  frequencies_grid.frequencies = (double *)malloc((frequencies_grid.N_frequencies) * sizeof(double));
  if (frequencies_grid.frequencies == NULL) {
    fprintf(stderr, "Error allocating memory for frequencies array\n");
    H5Gclose(freq_group);
    return frequencies_grid;
  }
  H5LTread_dataset_double(freq_group, "frequencies", frequencies_grid.frequencies);
  H5Gclose(freq_group);
  return frequencies_grid;
}

//============================================================================
// Empty structure builders
//============================================================================

/////////////////////////////////////////////////////////////
/// THDF_make_empty_field
/////////////////////////////////////////////////////////////
THDF_field_t
THDF_make_empty_field(void) {
  THDF_field_t output_field;
  output_field.index_i       = -1;
  output_field.index_j       = -1;
  output_field.index_incl    = -1;
  output_field.index_azimuth = -1;
  output_field.stokes_I      = NULL;
  output_field.stokes_QI     = NULL;
  output_field.stokes_UI     = NULL;
  output_field.stokes_VI     = NULL;

  output_field.norm_multiplier_I  = NULL;
  output_field.norm_multiplier_QI = NULL;
  output_field.norm_multiplier_UI = NULL;
  output_field.norm_multiplier_VI = NULL;

  return output_field;
}

/////////////////////////////////////////////////////////////
/// THDF_make_empty_angular_grid
/////////////////////////////////////////////////////////////
THDF_angular_grid_t
THDF_make_empty_angular_grid(void) {
  THDF_angular_grid_t angular_grid;
  angular_grid.N_directions         = -1;
  angular_grid.inclination_angles   = NULL;
  angular_grid.azimuthal_angles     = NULL;
  angular_grid.inclinations_indices = NULL;
  angular_grid.azimuthal_indices    = NULL;
  return angular_grid;
}

/////////////////////////////////////////////////////////////
/// THDF_make_empty_frequencies_grid
/////////////////////////////////////////////////////////////
THDF_frequencies_grid_t
THDF_make_empty_frequencies_grid(void) {
  THDF_frequencies_grid_t frequencies_grid;
  frequencies_grid.N_frequencies = 0;
  frequencies_grid.frequencies   = NULL;
  return frequencies_grid;
}

//============================================================================
// Clear / free helpers
//============================================================================

/////////////////////////////////////////////////////////////
/// THDF_clear_field
/////////////////////////////////////////////////////////////
void
THDF_clear_field(THDF_field_t *output_field) {

  if (!output_field)
    return;

  if (output_field->stokes_I)
    free(output_field->stokes_I);
  if (output_field->stokes_QI)
    free(output_field->stokes_QI);
  if (output_field->stokes_UI)
    free(output_field->stokes_UI);
  if (output_field->stokes_VI)
    free(output_field->stokes_VI);

  output_field->stokes_I           = NULL;
  output_field->stokes_QI          = NULL;
  output_field->stokes_UI          = NULL;
  output_field->stokes_VI          = NULL;
  output_field->norm_multiplier_I  = NULL;
  output_field->norm_multiplier_QI = NULL;
  output_field->norm_multiplier_UI = NULL;
  output_field->norm_multiplier_VI = NULL;
  output_field->index_i            = -1;
  output_field->index_j            = -1;
  output_field->index_incl         = -1;
  output_field->index_azimuth      = -1;
}

/////////////////////////////////////////////////////////////
/// THDF_clear_angular_grid
/////////////////////////////////////////////////////////////
void
THDF_clear_angular_grid(THDF_angular_grid_t *angular_grid) {
  if (!angular_grid)
    return;
  if (angular_grid->inclination_angles)
    free(angular_grid->inclination_angles);
  if (angular_grid->azimuthal_angles)
    free(angular_grid->azimuthal_angles);
  if (angular_grid->inclinations_indices)
    free(angular_grid->inclinations_indices);
  if (angular_grid->azimuthal_indices)
    free(angular_grid->azimuthal_indices);
  angular_grid->inclination_angles   = NULL;
  angular_grid->azimuthal_angles     = NULL;
  angular_grid->inclinations_indices = NULL;
  angular_grid->azimuthal_indices    = NULL;
  angular_grid->N_directions         = -1;
  angular_grid->N_inclination_angles = 0;
  angular_grid->N_azimuthal_angles   = 0;
}

/////////////////////////////////////////////////////////////
/// THDF_clear_frequencies_grid
/////////////////////////////////////////////////////////////
void
THDF_clear_frequencies_grid(THDF_frequencies_grid_t *frequencies_grid) {
  if (!frequencies_grid)
    return;
  if (frequencies_grid->frequencies)
    free(frequencies_grid->frequencies);
  frequencies_grid->frequencies   = NULL;
  frequencies_grid->N_frequencies = 0;
}
