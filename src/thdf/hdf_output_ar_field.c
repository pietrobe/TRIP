#include "hdf_output_ar_field.h"

// Close all HDF5 handles in the dataset handler if they are valid
#define CLOSE_OUTPUT_AR_FIELD(d)    \
  do {                              \
    if ((d).dataset_id_I >= 0)      \
      H5Dclose((d).dataset_id_I);   \
    if ((d).dataset_id_Q >= 0)      \
      H5Dclose((d).dataset_id_Q);   \
    if ((d).dataset_id_U >= 0)      \
      H5Dclose((d).dataset_id_U);   \
    if ((d).dataset_id_V >= 0)      \
      H5Dclose((d).dataset_id_V);   \
    if ((d).datatype_id >= 0)       \
      H5Tclose((d).datatype_id);    \
    if ((d).dataspace_id_I >= 0)    \
      H5Sclose((d).dataspace_id_I); \
    if ((d).dataspace_id_Q >= 0)    \
      H5Sclose((d).dataspace_id_Q); \
    if ((d).dataspace_id_U >= 0)    \
      H5Sclose((d).dataspace_id_U); \
    if ((d).dataspace_id_V >= 0)    \
      H5Sclose((d).dataspace_id_V); \
    if ((d).group_id >= 0)          \
      H5Gclose((d).group_id);       \
    (d).is_open = 0;                \
  } while (0)

/////////////////////////////////////////////////////////////
/// THDF_make_empty_ard_field_handler
/////////////////////////////////////////////////////////////
THDF_ard_field_handler
THDF_make_empty_ard_field_handler(void) {
  THDF_ard_field_handler output_dset;

  output_dset.file_id      = -1;
  output_dset.group_id     = -1;
  output_dset.dataset_id_I = -1;
  output_dset.dataset_id_Q = -1;
  output_dset.dataset_id_U = -1;
  output_dset.dataset_id_V = -1;
  // output_dset.dataset_id_mu     = -1;
  // output_dset.dataset_id_chi    = -1;
  output_dset.dataspace_id_I = -1;
  output_dset.dataspace_id_Q = -1;
  output_dset.dataspace_id_U = -1;
  output_dset.dataspace_id_V = -1;
  // output_dset.dataspace_mu      = -1;
  // output_dset.dataspace_chi     = -1;
  output_dset.datatype_id       = -1;
  output_dset.N_x               = 0;
  output_dset.N_y               = 0;
  output_dset.N_arbitrary_beams = 0;
  output_dset.N_frequencies     = 0;
  output_dset.is_open           = 0;

  return output_dset;
}  // END THDF_make_empty_ard_field_handler

/////////////////////////////////////////////////////////////
/// THDF_create_ard_field_handler_mpi
/////////////////////////////////////////////////////////////
THDF_ard_field_handler *
THDF_create_ard_field_handler_mpi(const hid_t  file_id,            //
                                  const double mu,                 //
                                  const double chi,                //
                                  const int    N_x,                //
                                  const int    N_y,                //
                                  const int    N_arbitrary_beams,  //
                                  const int    N_frequencies) {       //

  THDF_ard_field_handler *output_dset = malloc(sizeof(THDF_ard_field_handler));
  *output_dset                        = THDF_make_empty_ard_field_handler();

  output_dset->file_id           = file_id;
  output_dset->N_x               = N_x;
  output_dset->N_y               = N_y;
  output_dset->N_arbitrary_beams = N_arbitrary_beams;
  output_dset->N_frequencies     = N_frequencies;

  // Open or create group for arbitrary field output
  H5E_BEGIN_TRY { output_dset->group_id = H5Gopen2(file_id, "/emergent_beams", H5P_DEFAULT); }
  H5E_END_TRY;

  if (output_dset->group_id < 0) {
    // Group doesn't exist, create it
    output_dset->group_id = H5Gcreate2(file_id, "/emergent_beams", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (output_dset->group_id < 0) {
      fprintf(stderr, "Error creating emergent_Omega_beams group\n");
      output_dset->is_open = 0;
      free(output_dset);
      return NULL;
    }
  }

  // Define dataspace dimensions: [N_x, N_y, N_arbitrary_beams, N_frequencies]
  hsize_t dims[3] = {(hsize_t)N_x, (hsize_t)N_y, (hsize_t)N_frequencies};
  //   hsize_t dims_direction[4] = {(hsize_t)N_x, (hsize_t)N_y, (hsize_t)N_arbitrary_beams, (hsize_t)1};

  // Create datatype for Stokes parameters (copy predefined type so we can close it later)
  output_dset->datatype_id = H5Tcopy(THDF_get_hdf_float_datatype());
  if (output_dset->datatype_id < 0) {
    fprintf(stderr, "Error getting datatype for output_ar_field\n");
    CLOSE_OUTPUT_AR_FIELD((*output_dset));
    free(output_dset);
    return NULL;
  }

  // Create dataspace for each Stokes parameter dataset
  output_dset->dataspace_id_I = H5Screate_simple(3, dims, NULL);
  output_dset->dataspace_id_Q = H5Screate_simple(3, dims, NULL);
  output_dset->dataspace_id_U = H5Screate_simple(3, dims, NULL);
  output_dset->dataspace_id_V = H5Screate_simple(3, dims, NULL);

  // Create or open subgroups for each Stokes parameter
  hid_t group_I, group_Q, group_U, group_V;

  H5E_BEGIN_TRY { group_I = H5Gopen2(output_dset->group_id, "beams_I", H5P_DEFAULT); }
  H5E_END_TRY;
  if (group_I < 0) {
    group_I = H5Gcreate2(output_dset->group_id, "beams_I", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  H5E_BEGIN_TRY { group_Q = H5Gopen2(output_dset->group_id, "beams_QI_pc", H5P_DEFAULT); }
  H5E_END_TRY;
  if (group_Q < 0) {
    group_Q = H5Gcreate2(output_dset->group_id, "beams_QI_pc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  H5E_BEGIN_TRY { group_U = H5Gopen2(output_dset->group_id, "beams_UI_pc", H5P_DEFAULT); }
  H5E_END_TRY;
  if (group_U < 0) {
    group_U = H5Gcreate2(output_dset->group_id, "beams_UI_pc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  H5E_BEGIN_TRY { group_V = H5Gopen2(output_dset->group_id, "beams_VI_pc", H5P_DEFAULT); }
  H5E_END_TRY;
  if (group_V < 0) {
    group_V = H5Gcreate2(output_dset->group_id, "beams_VI_pc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  if (group_I < 0 || group_Q < 0 || group_U < 0 || group_V < 0) {

    fprintf(stderr, "Error creating Stokes parameter subgroups\n");

    if (group_I >= 0)
      H5Gclose(group_I);
    if (group_Q >= 0)
      H5Gclose(group_Q);
    if (group_U >= 0)
      H5Gclose(group_U);
    if (group_V >= 0)
      H5Gclose(group_V);

    CLOSE_OUTPUT_AR_FIELD((*output_dset));

    free(output_dset);

    return NULL;
  }  // END if error creating subgroups

  char dataset_name[256];
  snprintf(dataset_name, sizeof(dataset_name), "mu_%.3f_chi_%.3f", mu, chi);

  // Create datasets for Stokes parameters within their respective subgroups
  output_dset->dataset_id_I = H5Dcreate2(group_I, dataset_name, output_dset->datatype_id, output_dset->dataspace_id_I,
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  output_dset->dataset_id_Q = H5Dcreate2(group_Q, dataset_name, output_dset->datatype_id, output_dset->dataspace_id_Q,
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  output_dset->dataset_id_U = H5Dcreate2(group_U, dataset_name, output_dset->datatype_id, output_dset->dataspace_id_U,
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  output_dset->dataset_id_V = H5Dcreate2(group_V, dataset_name, output_dset->datatype_id, output_dset->dataspace_id_V,
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Close the subgroup handles (datasets remain accessible)
  H5Gclose(group_I);
  H5Gclose(group_Q);
  H5Gclose(group_U);
  H5Gclose(group_V);

  if (output_dset->dataset_id_I < 0 || output_dset->dataset_id_Q < 0 || output_dset->dataset_id_U < 0 ||
      output_dset->dataset_id_V < 0) {
    fprintf(stderr, "Error creating output_ar_field datasets\n");
    CLOSE_OUTPUT_AR_FIELD((*output_dset));
    free(output_dset);
    return NULL;
  }

  output_dset->is_open = 1;
  return output_dset;
}  // END THDF_create_ard_field_handler_mpi

/////////////////////////////////////////////////////////////
/// THDF_write_ard_field_to_hdf5
/////////////////////////////////////////////////////////////
int
THDF_write_ard_field_to_hdf5(THDF_ard_field_handler *output_dset,   //
                             THDF_ard_field_t       *output_field,  //
                             const hsize_t           start_i,       //
                             const hsize_t           start_j,       //
                             const hsize_t           count_i,       //
                             const hsize_t           count_j,       //
                             const hsize_t           count_frequencies) {

  // hsize_t mem_dims[3] = {count_i, count_j, count_frequencies};
  hid_t memspace = H5Screate_simple(3, (hsize_t[]){count_i, count_j, count_frequencies}, NULL);

  // Create separate dataspaces for each Stokes parameter
  hid_t dataspace_I = H5Dget_space(output_dset->dataset_id_I);
  hid_t dataspace_Q = H5Dget_space(output_dset->dataset_id_Q);
  hid_t dataspace_U = H5Dget_space(output_dset->dataset_id_U);
  hid_t dataspace_V = H5Dget_space(output_dset->dataset_id_V);

  hsize_t start[3] = {start_i, start_j, 0};
  hsize_t count[3] = {count_i, count_j, count_frequencies};

  H5Sselect_hyperslab(dataspace_I, H5S_SELECT_SET, start, NULL, count, NULL);
  herr_t status = H5Dwrite(output_dset->dataset_id_I, output_dset->datatype_id, memspace, dataspace_I, H5P_DEFAULT,
                           output_field->stokes_I);
  HDF5_CHECK_WRITE_OF(status, "Stokes I", memspace, -1);

  H5Sselect_hyperslab(dataspace_Q, H5S_SELECT_SET, start, NULL, count, NULL);
  status = H5Dwrite(output_dset->dataset_id_Q, output_dset->datatype_id, memspace, dataspace_Q, H5P_DEFAULT,
                    output_field->stokes_QI);
  HDF5_CHECK_WRITE_OF(status, "Stokes QI", memspace, -1);

  H5Sselect_hyperslab(dataspace_U, H5S_SELECT_SET, start, NULL, count, NULL);
  status = H5Dwrite(output_dset->dataset_id_U, output_dset->datatype_id, memspace, dataspace_U, H5P_DEFAULT,
                    output_field->stokes_UI);
  HDF5_CHECK_WRITE_OF(status, "Stokes UI", memspace, -1);

  H5Sselect_hyperslab(dataspace_V, H5S_SELECT_SET, start, NULL, count, NULL);
  status = H5Dwrite(output_dset->dataset_id_V, output_dset->datatype_id, memspace, dataspace_V, H5P_DEFAULT,
                    output_field->stokes_VI);
  HDF5_CHECK_WRITE_OF(status, "Stokes VI", memspace, -1);

  H5Sclose(memspace);
  H5Sclose(dataspace_I);
  H5Sclose(dataspace_Q);
  H5Sclose(dataspace_U);
  H5Sclose(dataspace_V);

  return 0;
}  // END THDF_write_ard_field_to_hdf5

/////////////////////////////////////////////////////////////
/// THDF_close_ard_field_handler_mpi
/////////////////////////////////////////////////////////////
void
THDF_close_ard_field_handler_mpi(THDF_ard_field_handler *output_dset) {
  if (!output_dset) {
    return;
  }

  if (output_dset->is_open) {
    CLOSE_OUTPUT_AR_FIELD(*output_dset);
  }

  output_dset->is_open = 0;

  free(output_dset);
}  // END THDF_close_ard_field_handler_mpi

/////////////////////////////////////////////////////////////
/// THDF_write_ard_directions
/////////////////////////////////////////////////////////////
int
THDF_write_ard_directions(hid_t                        file_id,           //
                          const THDF_ard_directions_t *ard_directions) {  //

  hid_t dir_group = H5Gcreate2(file_id, "/beams_directions", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dir_group < 0) {
    fprintf(stderr, "Error creating output arbitrary directions group\n");
    return -1;
  }

  H5LTmake_dataset_double(dir_group, "mu", 1, (hsize_t[]){ard_directions->N_directions}, ard_directions->mu);
  H5LTmake_dataset_double(dir_group, "chi", 1, (hsize_t[]){ard_directions->N_directions}, ard_directions->chi);
  H5LTmake_dataset_int(dir_group, "N_directions", 1, (hsize_t[]){1}, &ard_directions->N_directions);

  H5Gclose(dir_group);
  return 0;
}  // END THDF_write_ard_directions

/////////////////////////////////////////////////////////////
/// THDF_write_ard_frequencies
/////////////////////////////////////////////////////////////
int
THDF_write_ard_frequencies(hid_t                         file_id,      //
                           const THDF_ard_frequencies_t *ard_freqs) {  //

  hid_t freq_group = H5Gcreate2(file_id, "/frequencies_grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (freq_group < 0) {
    fprintf(stderr, "Error creating output arbitrary frequencies group\n");
    return -1;
  }

  H5LTmake_dataset_double(freq_group, "frequencies", 1, (hsize_t[]){ard_freqs->N_frequencies}, ard_freqs->frequencies);
  H5LTmake_dataset_int(freq_group, "N_frequencies", 1, (hsize_t[]){1}, &ard_freqs->N_frequencies);

  H5Gclose(freq_group);
  return 0;
}  // END THDF_write_ard_frequencies
