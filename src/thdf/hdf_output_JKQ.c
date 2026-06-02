#include <stdlib.h>

#include "hdf_output_JKQ.h"

#define HDF5_CHECK_STATUS(_status_, _msg_, _handler_)                              \
  if (_status_ < 0) {                                                              \
    char _buffer_[2048];                                                           \
    snprintf(_buffer_, 2048, "%s, file %s, line %d\n", _msg_, __FILE__, __LINE__); \
    fprintf(stderr, "%s", _buffer_);                                               \
    free(_handler_);                                                               \
    return NULL;                                                                   \
  }

//////////////////////////////////////////////////////////////
/// THDF_create_JKQ_field_handler_mpi
//////////////////////////////////////////////////////////////
THDF_JKQ_handler_t *                                      //
THDF_create_JKQ_field_handler_mpi(hid_t file,             //
                                  int   N_x,              //
                                  int   N_y,              //
                                  int   N_z,              //
                                  int   N_JKQ_real,       //
                                  int   N_JKQ_imag,       //
                                  int   N_frequencies) {  //

  THDF_JKQ_handler_t *handler = (THDF_JKQ_handler_t *)malloc(sizeof(THDF_JKQ_handler_t));
  handler->file_id            = file;
  handler->N_x                = N_x;
  handler->N_y                = N_y;
  handler->N_z                = N_z;
  handler->N_JKQ_real         = N_JKQ_real;
  handler->N_JKQ_imag         = N_JKQ_imag;
  handler->N_frequencies      = N_frequencies;

  handler->datatype_id           = H5T_NATIVE_FLOAT;
  handler->datatype_norm_mult_id = H5T_NATIVE_DOUBLE;

  handler->is_open = 1;

  hsize_t dims_re[5] = {N_x, N_y, N_z, N_JKQ_real, N_frequencies};
  hsize_t dims_im[5] = {N_x, N_y, N_z, N_JKQ_imag, N_frequencies};

  hsize_t dims_norm_re[4] = {N_x, N_y, N_z, N_JKQ_real};
  hsize_t dims_norm_im[4] = {N_x, N_y, N_z, N_JKQ_imag};

  // Create group for JKQ datasets (or open if it already exists)
  // Suppress errors temporarily to allow "group already exists" errors
  H5E_auto_t old_func;
  void      *old_client_data;
  H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  handler->group_id = H5Gcreate2(file, "/JKQ", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Restore error handling
  H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

  if (handler->group_id < 0) {
    // Group already exists, open it
    handler->group_id = H5Gopen2(file,          //
                                 "/JKQ",        //
                                 H5P_DEFAULT);  //
    HDF5_CHECK_STATUS(handler->group_id, "Error opening existing JKQ group in HDF5 file\n", handler);

    // Open existing datasets
    handler->dataset_id_JKQ_re = H5Dopen2(handler->group_id,  //
                                          "JKQ_re",           //
                                          H5P_DEFAULT);       //
    HDF5_CHECK_STATUS(handler->dataset_id_JKQ_re, "Error opening JKQ_re dataset in HDF5 file\n", handler);

    handler->dataset_id_JKQ_im = H5Dopen2(handler->group_id,  //
                                          "JKQ_im",           //
                                          H5P_DEFAULT);       //
    HDF5_CHECK_STATUS(handler->dataset_id_JKQ_im, "Error opening JKQ_im dataset in HDF5 file\n", handler);

    handler->dataset_JKQ_norm_mult_re = H5Dopen2(handler->group_id,  //
                                                 "JKQ_norm_re",      //
                                                 H5P_DEFAULT);       //
    HDF5_CHECK_STATUS(handler->dataset_JKQ_norm_mult_re, "Error opening JKQ_norm_re dataset in HDF5 file\n", handler);

    handler->dataset_JKQ_norm_mult_im = H5Dopen2(handler->group_id,  //
                                                 "JKQ_norm_im",      //
                                                 H5P_DEFAULT);       //
    HDF5_CHECK_STATUS(handler->dataset_JKQ_norm_mult_im, "Error opening JKQ_norm_im dataset in HDF5 file\n", handler);

    // Retrieve dataspaces
    handler->dataspace_id_JKQ_re = H5Dget_space(handler->dataset_id_JKQ_re);  //
    HDF5_CHECK_STATUS(handler->dataspace_id_JKQ_re, "Error getting dataspace for JKQ_re dataset\n", handler);

    handler->dataspace_id_JKQ_im = H5Dget_space(handler->dataset_id_JKQ_im);  //
    HDF5_CHECK_STATUS(handler->dataspace_id_JKQ_im, "Error getting dataspace for JKQ_im dataset\n", handler);

    handler->dataspace_JKQ_norm_mult_re = H5Dget_space(handler->dataset_JKQ_norm_mult_re);  //
    HDF5_CHECK_STATUS(handler->dataspace_JKQ_norm_mult_re, "Error getting dataspace for JKQ_norm_re dataset\n", handler);

    handler->dataspace_JKQ_norm_mult_im = H5Dget_space(handler->dataset_JKQ_norm_mult_im);  //
    HDF5_CHECK_STATUS(handler->dataspace_JKQ_norm_mult_im, "Error getting dataspace for JKQ_norm_im dataset\n", handler);
  } else {
    // Group was created, create the datasets
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);

    // Create datasets for JKQ real part
    handler->dataspace_id_JKQ_re = H5Screate_simple(5, dims_re, NULL);
    handler->dataset_id_JKQ_re   = H5Dcreate2(handler->group_id,             //
                                              "JKQ_re",                      //
                                              handler->datatype_id,          //
                                              handler->dataspace_id_JKQ_re,  //
                                              H5P_DEFAULT,                   //
                                              dcpl_id,                       //
                                              H5P_DEFAULT);                  //
    HDF5_CHECK_STATUS(handler->dataset_id_JKQ_re, "Error creating JKQ_re dataset in HDF5 file\n", handler);

    // Create datasets for JKQ imaginary part
    handler->dataspace_id_JKQ_im = H5Screate_simple(5, dims_im, NULL);
    handler->dataset_id_JKQ_im   = H5Dcreate2(handler->group_id,             //
                                              "JKQ_im",                      //
                                              handler->datatype_id,          //
                                              handler->dataspace_id_JKQ_im,  //
                                              H5P_DEFAULT,                   //
                                              dcpl_id,                       //
                                              H5P_DEFAULT);                  //
    HDF5_CHECK_STATUS(handler->dataset_id_JKQ_im, "Error creating JKQ_im dataset in HDF5 file\n", handler);

    // Create datasets for normalization multipliers (real part)
    handler->dataspace_JKQ_norm_mult_re = H5Screate_simple(4, dims_norm_re, NULL);
    handler->dataset_JKQ_norm_mult_re   = H5Dcreate2(handler->group_id,                    //
                                                     "JKQ_norm_re",                        //
                                                     handler->datatype_norm_mult_id,       //
                                                     handler->dataspace_JKQ_norm_mult_re,  //
                                                     H5P_DEFAULT,                          //
                                                     dcpl_id,                              //
                                                     H5P_DEFAULT);                         //
    HDF5_CHECK_STATUS(handler->dataset_JKQ_norm_mult_re, "Error creating JKQ_norm_re dataset in HDF5 file\n", handler);

    // Create datasets for normalization multipliers (imaginary part)
    handler->dataspace_JKQ_norm_mult_im = H5Screate_simple(4, dims_norm_im, NULL);
    handler->dataset_JKQ_norm_mult_im   = H5Dcreate2(handler->group_id,                    //
                                                     "JKQ_norm_im",                        //
                                                     handler->datatype_norm_mult_id,       //
                                                     handler->dataspace_JKQ_norm_mult_im,  //
                                                     H5P_DEFAULT,                          //
                                                     dcpl_id,                              //
                                                     H5P_DEFAULT);                         //
    HDF5_CHECK_STATUS(handler->dataset_JKQ_norm_mult_im, "Error creating JKQ_norm_im dataset in HDF5 file\n", handler);

    H5Pclose(dcpl_id);
  }

  return handler;
}

//////////////////////////////////////////////////////////////
/// THDF_close_JKQ_field_handler_mpi
//////////////////////////////////////////////////////////////
void                                                             //
THDF_close_JKQ_field_handler_mpi(THDF_JKQ_handler_t *handler) {  //

  if (handler == NULL) {
    return;
  }

  H5Dclose(handler->dataset_id_JKQ_re);
  H5Sclose(handler->dataspace_id_JKQ_re);

  H5Dclose(handler->dataset_id_JKQ_im);
  H5Sclose(handler->dataspace_id_JKQ_im);

  H5Dclose(handler->dataset_JKQ_norm_mult_re);
  H5Sclose(handler->dataspace_JKQ_norm_mult_re);

  H5Dclose(handler->dataset_JKQ_norm_mult_im);
  H5Sclose(handler->dataspace_JKQ_norm_mult_im);

  H5Gclose(handler->group_id);

  free(handler);
}

#define HDF5_CHECK_KQ_STATUS(_status_, _msg_, _ret_)                                 \
  if (_status_ < 0) {                                                                \
    char __buffer__[2048];                                                           \
    snprintf(__buffer__, 2048, "%s, file %s, line %d\n", _msg_, __FILE__, __LINE__); \
    fprintf(stderr, "%s", __buffer__);                                               \
    return _ret_;                                                                    \
  }

//////////////////////////////////////////////////////////////
/// THDF_write_KQ_table_to_hdf5
//////////////////////////////////////////////////////////////
int
THDF_write_KQ_table_to_hdf5(hid_t            file,        //
                            THDF_KQ_table_t *KQ_table) {  //

  // Create group for KQ matrix
  hid_t group_id = H5Gcreate2(file, "/KQ_table", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_CHECK_KQ_STATUS(group_id, "Error creating KQ_table group in HDF5 file\n", -1);

  // Create dataset for KQ matrix
  hsize_t dims_KQ[2]      = {KQ_table->KQ_size, 2};  // 2 columns: K and Q
  hid_t   dataspace_id_KQ = H5Screate_simple(2, dims_KQ, NULL);
  hid_t   dataset_id_KQ = H5Dcreate2(group_id,  //
                                     "KQ_table", H5T_NATIVE_INT, dataspace_id_KQ, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_CHECK_KQ_STATUS(dataset_id_KQ, "Error creating KQ_table dataset in HDF5 file\n", -1);

  // Write KQ matrix data
  herr_t status = H5Dwrite(dataset_id_KQ, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, KQ_table->KQ_table);
  HDF5_CHECK_KQ_STATUS(status, "Error writing KQ_table data to HDF5 file\n", -1);

  // Close resources
  H5Dclose(dataset_id_KQ);
  H5Sclose(dataspace_id_KQ);

  hsize_t dims_comp_KQ_real_comm[2] = {KQ_table->KQ_compressed_size_real, 2};
  hsize_t dims_comp_KQ_imag_comm[2] = {KQ_table->KQ_compressed_size_imag, 2};

  // Create dataset for compressed KQ real part
  hid_t dataspace_id_comp_KQ_real = H5Screate_simple(2, dims_comp_KQ_real_comm, NULL);
  hid_t dataset_id_comp_KQ_real   = H5Dcreate2(group_id, "KQ_compressed_real", H5T_NATIVE_INT, dataspace_id_comp_KQ_real,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_CHECK_KQ_STATUS(dataset_id_comp_KQ_real, "Error creating KQ_compressed_real dataset in HDF5 file\n", -1);

  // Write compressed KQ real part data
  status = H5Dwrite(dataset_id_comp_KQ_real,                                                       //
                    H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, KQ_table->KQ_compressed_real);  //
  HDF5_CHECK_KQ_STATUS(status, "Error writing KQ_compressed_real data to HDF5 file\n", -1);

  // Create dataset for compressed KQ imaginary part
  hid_t dataspace_id_comp_KQ_imag = H5Screate_simple(2, dims_comp_KQ_imag_comm, NULL);
  hid_t dataset_id_comp_KQ_imag   = H5Dcreate2(group_id, "KQ_compressed_imag", H5T_NATIVE_INT, dataspace_id_comp_KQ_imag,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_CHECK_KQ_STATUS(dataset_id_comp_KQ_imag, "Error creating KQ_compressed_imag dataset in HDF5 file\n", -1);

  // Write compressed KQ imaginary part data
  status = H5Dwrite(dataset_id_comp_KQ_imag,                                                       //
                    H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, KQ_table->KQ_compressed_imag);  //
  HDF5_CHECK_KQ_STATUS(status, "Error writing KQ_compressed_imag data to HDF5 file\n", -1);

  status = H5Dwrite(dataset_id_comp_KQ_imag,                                                       //
                    H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, KQ_table->KQ_compressed_imag);  //
  HDF5_CHECK_KQ_STATUS(status, "Error writing KQ_compressed_imag data to HDF5 file\n", -1);

  // Close resources
  H5Dclose(dataset_id_comp_KQ_real);
  H5Sclose(dataspace_id_comp_KQ_real);

  H5Dclose(dataset_id_comp_KQ_imag);
  H5Sclose(dataspace_id_comp_KQ_imag);

  H5Gclose(group_id);

  return 0;
}

#define HDF5_CHECK_JKQ_STATUS(_status_, _msg_, _ret_)                                \
  if (_status_ < 0) {                                                                \
    char __buffer__[2048];                                                           \
    snprintf(__buffer__, 2048, "%s, file %s, line %d\n", _msg_, __FILE__, __LINE__); \
    fprintf(stderr, "%s", __buffer__);                                               \
    return _ret_;                                                                    \
  }

//////////////////////////////////////////////////////////////
/// THDF_write_JKQ_field_to_hdf5
//////////////////////////////////////////////////////////////
int
THDF_write_JKQ_field_to_hdf5(THDF_JKQ_handler_t *handler,    //
                             THDF_JKQ_field_t   *JKQ_field,  //
                             hsize_t             start_i,    //
                             hsize_t             start_j,    //
                             hsize_t             start_k,    //
                             hsize_t             count_i,    //
                             hsize_t             count_j,    //
                             hsize_t             count_k) {  //

  // Define hyperslab in the file
  hsize_t file_start[5]    = {start_i, start_j, start_k, 0, 0};
  hsize_t file_count_re[5] = {count_i, count_j, count_k, handler->N_JKQ_real, handler->N_frequencies};
  hsize_t file_count_im[5] = {count_i, count_j, count_k, handler->N_JKQ_imag, handler->N_frequencies};

  // Create transfer property list for collective I/O
  hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);

  herr_t status;
  // Write JKQ real part
  status = H5Sselect_hyperslab(handler->dataspace_id_JKQ_re, H5S_SELECT_SET, file_start, NULL, file_count_re, NULL);
  HDF5_CHECK_JKQ_STATUS(status, "Error selecting hyperslab for JKQ_re dataset\n", -1);

  hid_t memspace_id_re = H5Screate_simple(5, file_count_re, NULL);
  status               = H5Dwrite(handler->dataset_id_JKQ_re,    //
                                  handler->datatype_id,          //
                                  memspace_id_re,                //
                                  handler->dataspace_id_JKQ_re,  //
                                  dxpl_id,                       //
                                  JKQ_field->JKQ_re);            //
  HDF5_CHECK_JKQ_STATUS(status, "Error writing JKQ_re data to HDF5 file\n", -1);

  H5Sclose(memspace_id_re);

  // Write JKQ imaginary part
  status = H5Sselect_hyperslab(handler->dataspace_id_JKQ_im, H5S_SELECT_SET, file_start, NULL, file_count_im, NULL);
  HDF5_CHECK_JKQ_STATUS(status, "Error selecting hyperslab for JKQ_im dataset\n", -1);

  hid_t memspace_id_im = H5Screate_simple(5, file_count_im, NULL);
  status               = H5Dwrite(handler->dataset_id_JKQ_im,    //
                                  handler->datatype_id,          //
                                  memspace_id_im,                //
                                  handler->dataspace_id_JKQ_im,  //
                                  dxpl_id,                       //
                                  JKQ_field->JKQ_im);            //
  HDF5_CHECK_JKQ_STATUS(status, "Error writing JKQ_im data to HDF5 file\n", -1);

  H5Sclose(memspace_id_im);

  // Write normalization multipliers (real part)
  hsize_t file_start_norm_re[4] = {start_i, start_j, start_k, 0};
  hsize_t file_count_norm_re[4] = {count_i, count_j, count_k, handler->N_JKQ_real};

  status = H5Sselect_hyperslab(handler->dataspace_JKQ_norm_mult_re, H5S_SELECT_SET, file_start_norm_re, NULL,
                               file_count_norm_re, NULL);
  HDF5_CHECK_JKQ_STATUS(status, "Error selecting hyperslab for JKQ_norm_re dataset\n", -1);

  hid_t memspace_id_norm_re = H5Screate_simple(4, file_count_norm_re, NULL);
  status                    = H5Dwrite(handler->dataset_JKQ_norm_mult_re,    //
                                       handler->datatype_norm_mult_id,       //
                                       memspace_id_norm_re,                  //
                                       handler->dataspace_JKQ_norm_mult_re,  //
                                       dxpl_id,                              //
                                       JKQ_field->norm_multiplier_JKQ_re);   //
  HDF5_CHECK_JKQ_STATUS(status, "Error writing JKQ_norm_re data to HDF5 file\n", -1);

  H5Sclose(memspace_id_norm_re);

  // Write normalization multipliers (imaginary part)
  hsize_t file_start_norm_im[4] = {start_i, start_j, start_k, 0};
  hsize_t file_count_norm_im[4] = {count_i, count_j, count_k, handler->N_JKQ_imag};

  status = H5Sselect_hyperslab(handler->dataspace_JKQ_norm_mult_im, H5S_SELECT_SET, file_start_norm_im, NULL,
                               file_count_norm_im, NULL);
  HDF5_CHECK_JKQ_STATUS(status, "Error selecting hyperslab for JKQ_norm_im dataset\n", -1);

  hid_t memspace_id_norm_im = H5Screate_simple(4, file_count_norm_im, NULL);
  status                    = H5Dwrite(handler->dataset_JKQ_norm_mult_im,    //
                                       handler->datatype_norm_mult_id,       //
                                       memspace_id_norm_im,                  //
                                       handler->dataspace_JKQ_norm_mult_im,  //
                                       dxpl_id,                              //
                                       JKQ_field->norm_multiplier_JKQ_im);
  HDF5_CHECK_JKQ_STATUS(status, "Error writing JKQ_norm_im data to HDF5 file\n", -1);

  H5Sclose(memspace_id_norm_im);

  H5Pclose(dxpl_id);

  return 0;
}