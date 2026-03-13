#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf_output_3D_field.h"

/////////////////////////////////////////////////////
// write_stokes_hyperslab
/////////////////////////////////////////////////////
static int
write_stokes_hyperslab(hid_t          dataset_id,   //
                       hid_t          datatype_id,  //
                       hid_t          memspace,     //
                       const hsize_t *start,        //
                       const hsize_t *count,        //
                       const void    *data,         //
                       const char    *stokes_name) {   //

  hid_t dataspace = H5Dget_space(dataset_id);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, NULL);
  herr_t status = H5Dwrite(dataset_id, datatype_id, memspace, dataspace, H5P_DEFAULT, data);
  HDF5_CHECK_WRITE_OF(status, stokes_name, memspace, -1);
  H5Sclose(dataspace);

  return (int)status;
}

/////////////////////////////////////////////////////
// get_input_index
/////////////////////////////////////////////////////
static inline int
get_input_index(const int i,                      //
                const int j,                      //
                const int k,                      //
                const int incl,                   //
                const int az,                     //
                const int f,                      //
                const int stokes_offset,          //
                const int N0_in_local_x,          //
                const int N0_in_local_y,          //
                const int N0_in_local_z,          //
                const int stride_in_x,            //
                const int stride_in_y,            //
                const int stride_in_z,            //
                const int stride_in_incl,         //
                const int stride_in_azimuth,      //
                const int stride_in_frequencies,  //
                const int stride_in_stokes) {     //

  return (i + N0_in_local_x) * stride_in_x + (j + N0_in_local_y) * stride_in_y + (k + N0_in_local_z) * stride_in_z +
         incl * stride_in_incl + az * stride_in_azimuth + f * stride_in_frequencies + stokes_offset * stride_in_stokes;
}

/////////////////////////////////////////////////////
// THDF_create_3D_field_handler_mpi
/////////////////////////////////////////////////////
THDF_3D_field_handler_t *
THDF_create_3D_field_handler_mpi(hid_t      file,               //
                                 const bool normalized_output,  //
                                 const int  N_x,                //
                                 const int  N_y,                //
                                 const int  N_z,                //
                                 const int  N_incl,             //
                                 const int  N_azimuth,          //
                                 const int  N_frequencies) {     //

  THDF_3D_field_handler_t *output_dset = (THDF_3D_field_handler_t *)malloc(sizeof(THDF_3D_field_handler_t));

  if (!output_dset) {
    fprintf(stderr, "Error allocating memory for output field dataset struct\n");
    return NULL;
  }

  output_dset->file_id                      = -1;
  output_dset->group_id                     = -1;
  output_dset->dataset_id_I                 = -1;
  output_dset->dataset_id_Q                 = -1;
  output_dset->dataset_id_U                 = -1;
  output_dset->dataset_id_V                 = -1;
  output_dset->dataspace_id_I               = -1;
  output_dset->dataspace_id_Q               = -1;
  output_dset->dataspace_id_U               = -1;
  output_dset->dataspace_id_V               = -1;
  output_dset->dataset_norm_multiplier_I    = -1;
  output_dset->dataset_norm_multiplier_QI   = -1;
  output_dset->dataset_norm_multiplier_UI   = -1;
  output_dset->dataset_norm_multiplier_VI   = -1;
  output_dset->dataspace_norm_multiplier_I  = -1;
  output_dset->dataspace_norm_multiplier_QI = -1;
  output_dset->dataspace_norm_multiplier_UI = -1;
  output_dset->dataspace_norm_multiplier_VI = -1;
  output_dset->datatype_id                  = -1;
  output_dset->datatype_f32_id              = -1;
  output_dset->is_open                      = false;

  output_dset->file_id           = file;
  output_dset->N_x               = N_x;
  output_dset->N_y               = N_y;
  output_dset->N_z               = N_z;
  output_dset->N_incl            = N_incl;
  output_dset->N_azimuth         = N_azimuth;
  output_dset->N_frequencies     = N_frequencies;
  output_dset->normalized_output = normalized_output;

  hsize_t dims_field[6] = {N_x, N_y, N_z, N_incl, N_azimuth, N_frequencies};
  hsize_t dims_norm[5]  = {N_x, N_y, N_z, N_incl, N_azimuth};

  // Check if output_field group already exists, if so, open it
  bool group_existed = false;
  H5E_BEGIN_TRY { output_dset->group_id = H5Gopen2(file, "/radiation_field", H5P_DEFAULT); }
  H5E_END_TRY;

  // If group doesn't exist, create it (collective operation)
  if (output_dset->group_id < 0) {
    output_dset->group_id = H5Gcreate2(file, "/radiation_field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (output_dset->group_id < 0) {
      fprintf(stderr, "Error creating radiation_field group with MPI\n");
      free(output_dset);
      return NULL;
    }
  } else {
    group_existed = true;
  }

  // Create the datatype for float 32 bits
  output_dset->datatype_f32_id = H5Tcopy(THDF_get_hdf_float32_datatype());

  if (group_existed) {
    // Datasets already exist — open them instead of creating
    output_dset->dataset_id_I = H5Dopen2(output_dset->group_id, "I", H5P_DEFAULT);
    output_dset->dataset_id_Q = H5Dopen2(output_dset->group_id, "QI_pc", H5P_DEFAULT);
    output_dset->dataset_id_U = H5Dopen2(output_dset->group_id, "UI_pc", H5P_DEFAULT);
    output_dset->dataset_id_V = H5Dopen2(output_dset->group_id, "VI_pc", H5P_DEFAULT);

    output_dset->dataspace_id_I = H5Dget_space(output_dset->dataset_id_I);
    output_dset->dataspace_id_Q = H5Dget_space(output_dset->dataset_id_Q);
    output_dset->dataspace_id_U = H5Dget_space(output_dset->dataset_id_U);
    output_dset->dataspace_id_V = H5Dget_space(output_dset->dataset_id_V);
  } else {
    // Group is new — create datasets
    output_dset->dataspace_id_I = H5Screate_simple(6, dims_field, NULL);
    output_dset->dataset_id_I   = H5Dcreate2(output_dset->group_id, "I", output_dset->datatype_f32_id,
                                             output_dset->dataspace_id_I, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    output_dset->dataspace_id_Q = H5Screate_simple(6, dims_field, NULL);
    output_dset->dataset_id_Q   = H5Dcreate2(output_dset->group_id, "QI_pc", output_dset->datatype_f32_id,
                                             output_dset->dataspace_id_Q, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    output_dset->dataspace_id_U = H5Screate_simple(6, dims_field, NULL);
    output_dset->dataset_id_U   = H5Dcreate2(output_dset->group_id, "UI_pc", output_dset->datatype_f32_id,
                                             output_dset->dataspace_id_U, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    output_dset->dataspace_id_V = H5Screate_simple(6, dims_field, NULL);
    output_dset->dataset_id_V   = H5Dcreate2(output_dset->group_id, "VI_pc", output_dset->datatype_f32_id,
                                             output_dset->dataspace_id_V, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  if (output_dset->dataset_id_I < 0 || output_dset->dataset_id_Q < 0 || output_dset->dataset_id_U < 0 ||
      output_dset->dataset_id_V < 0) {
    fprintf(stderr, "Error creating output field dataset with MPI\n");
    H5Tclose(output_dset->datatype_f32_id);
    H5Sclose(output_dset->dataspace_id_I);
    H5Sclose(output_dset->dataspace_id_Q);
    H5Sclose(output_dset->dataspace_id_U);
    H5Sclose(output_dset->dataspace_id_V);
    H5Dclose(output_dset->dataset_id_I);
    H5Dclose(output_dset->dataset_id_Q);
    H5Dclose(output_dset->dataset_id_U);
    H5Dclose(output_dset->dataset_id_V);

    H5Gclose(output_dset->group_id);
    free(output_dset);
    return NULL;
  }

  if (!normalized_output) {
    // If not normalized output, we don't need to create the normalization multiplier datasets
    output_dset->dataset_norm_multiplier_I  = -1;
    output_dset->dataset_norm_multiplier_QI = -1;
    output_dset->dataset_norm_multiplier_UI = -1;
    output_dset->dataset_norm_multiplier_VI = -1;

    output_dset->dataspace_norm_multiplier_I  = -1;
    output_dset->dataspace_norm_multiplier_QI = -1;
    output_dset->dataspace_norm_multiplier_UI = -1;
    output_dset->dataspace_norm_multiplier_VI = -1;

    output_dset->is_open = true;
    return output_dset;
  }

  output_dset->datatype_id = H5Tcopy(THDF_get_hdf_float_datatype());

  if (group_existed) {
    // Normalization datasets already exist — open them
    output_dset->dataset_norm_multiplier_I  = H5Dopen2(output_dset->group_id, "norm_multiplier_I", H5P_DEFAULT);
    output_dset->dataset_norm_multiplier_QI = H5Dopen2(output_dset->group_id, "norm_multiplier_QI_pc", H5P_DEFAULT);
    output_dset->dataset_norm_multiplier_UI = H5Dopen2(output_dset->group_id, "norm_multiplier_UI_pc", H5P_DEFAULT);
    output_dset->dataset_norm_multiplier_VI = H5Dopen2(output_dset->group_id, "norm_multiplier_VI_pc", H5P_DEFAULT);

    output_dset->dataspace_norm_multiplier_I  = H5Dget_space(output_dset->dataset_norm_multiplier_I);
    output_dset->dataspace_norm_multiplier_QI = H5Dget_space(output_dset->dataset_norm_multiplier_QI);
    output_dset->dataspace_norm_multiplier_UI = H5Dget_space(output_dset->dataset_norm_multiplier_UI);
    output_dset->dataspace_norm_multiplier_VI = H5Dget_space(output_dset->dataset_norm_multiplier_VI);
  } else {
    // Group is new — create normalization datasets
    output_dset->dataspace_norm_multiplier_I = H5Screate_simple(5, dims_norm, NULL);
    output_dset->dataset_norm_multiplier_I =
        H5Dcreate2(output_dset->group_id, "norm_multiplier_I", output_dset->datatype_id,
                   output_dset->dataspace_norm_multiplier_I, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    output_dset->dataspace_norm_multiplier_QI = H5Screate_simple(5, dims_norm, NULL);
    output_dset->dataset_norm_multiplier_QI =
        H5Dcreate2(output_dset->group_id, "norm_multiplier_QI_pc", output_dset->datatype_id,
                   output_dset->dataspace_norm_multiplier_QI, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    output_dset->dataspace_norm_multiplier_UI = H5Screate_simple(5, dims_norm, NULL);
    output_dset->dataset_norm_multiplier_UI =
        H5Dcreate2(output_dset->group_id, "norm_multiplier_UI_pc", output_dset->datatype_id,
                   output_dset->dataspace_norm_multiplier_UI, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    output_dset->dataspace_norm_multiplier_VI = H5Screate_simple(5, dims_norm, NULL);
    output_dset->dataset_norm_multiplier_VI =
        H5Dcreate2(output_dset->group_id, "norm_multiplier_VI_pc", output_dset->datatype_id,
                   output_dset->dataspace_norm_multiplier_VI, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  if (output_dset->dataset_norm_multiplier_I < 0 || output_dset->dataset_norm_multiplier_QI < 0 ||
      output_dset->dataset_norm_multiplier_UI < 0 || output_dset->dataset_norm_multiplier_VI < 0) {
    fprintf(stderr, "Error creating normalization multiplier datasets with MPI\n");
    if (output_dset->dataset_norm_multiplier_I >= 0) {
      H5Dclose(output_dset->dataset_norm_multiplier_I);
    }
    if (output_dset->dataset_norm_multiplier_QI >= 0) {
      H5Dclose(output_dset->dataset_norm_multiplier_QI);
    }
    if (output_dset->dataset_norm_multiplier_UI >= 0) {
      H5Dclose(output_dset->dataset_norm_multiplier_UI);
    }
    if (output_dset->dataset_norm_multiplier_VI >= 0) {
      H5Dclose(output_dset->dataset_norm_multiplier_VI);
    }
    if (output_dset->dataspace_norm_multiplier_I >= 0) {
      H5Sclose(output_dset->dataspace_norm_multiplier_I);
    }
    if (output_dset->dataspace_norm_multiplier_QI >= 0) {
      H5Sclose(output_dset->dataspace_norm_multiplier_QI);
    }
    if (output_dset->dataspace_norm_multiplier_UI >= 0) {
      H5Sclose(output_dset->dataspace_norm_multiplier_UI);
    }
    if (output_dset->dataspace_norm_multiplier_VI >= 0) {
      H5Sclose(output_dset->dataspace_norm_multiplier_VI);
    }
    if (output_dset->dataset_id_I >= 0) {
      H5Dclose(output_dset->dataset_id_I);
    }
    if (output_dset->dataset_id_Q >= 0) {
      H5Dclose(output_dset->dataset_id_Q);
    }
    if (output_dset->dataset_id_U >= 0) {
      H5Dclose(output_dset->dataset_id_U);
    }
    if (output_dset->dataset_id_V >= 0) {
      H5Dclose(output_dset->dataset_id_V);
    }
    if (output_dset->datatype_id >= 0) {
      H5Tclose(output_dset->datatype_id);
    }
    if (output_dset->datatype_f32_id >= 0) {
      H5Tclose(output_dset->datatype_f32_id);
    }
    if (output_dset->dataspace_id_I >= 0) {
      H5Sclose(output_dset->dataspace_id_I);
    }
    if (output_dset->dataspace_id_Q >= 0) {
      H5Sclose(output_dset->dataspace_id_Q);
    }
    if (output_dset->dataspace_id_U >= 0) {
      H5Sclose(output_dset->dataspace_id_U);
    }
    if (output_dset->dataspace_id_V >= 0) {
      H5Sclose(output_dset->dataspace_id_V);
    }
    if (output_dset->group_id >= 0) {
      H5Gclose(output_dset->group_id);
    }
    free(output_dset);
    return NULL;
  }

  output_dset->is_open = true;
  return output_dset;
}

/////////////////////////////////////////////////////
// THDF_close_3D_field_handler_mpi
/////////////////////////////////////////////////////
void
THDF_close_3D_field_handler_mpi(THDF_3D_field_handler_t *output_dset) {
  if (!output_dset) {
    return;
  }

  if (output_dset->is_open) {
    if (output_dset->dataset_id_I >= 0) {
      H5Dclose(output_dset->dataset_id_I);
    }
    if (output_dset->dataset_id_Q >= 0) {
      H5Dclose(output_dset->dataset_id_Q);
    }
    if (output_dset->dataset_id_U >= 0) {
      H5Dclose(output_dset->dataset_id_U);
    }
    if (output_dset->dataset_id_V >= 0) {
      H5Dclose(output_dset->dataset_id_V);
    }
    if (output_dset->dataspace_id_I >= 0) {
      H5Sclose(output_dset->dataspace_id_I);
    }
    if (output_dset->dataspace_id_Q >= 0) {
      H5Sclose(output_dset->dataspace_id_Q);
    }
    if (output_dset->dataspace_id_U >= 0) {
      H5Sclose(output_dset->dataspace_id_U);
    }
    if (output_dset->dataspace_id_V >= 0) {
      H5Sclose(output_dset->dataspace_id_V);
    }

    if (output_dset->dataset_norm_multiplier_I >= 0) {
      H5Dclose(output_dset->dataset_norm_multiplier_I);
    }
    if (output_dset->dataset_norm_multiplier_QI >= 0) {
      H5Dclose(output_dset->dataset_norm_multiplier_QI);
    }
    if (output_dset->dataset_norm_multiplier_UI >= 0) {
      H5Dclose(output_dset->dataset_norm_multiplier_UI);
    }
    if (output_dset->dataset_norm_multiplier_VI >= 0) {
      H5Dclose(output_dset->dataset_norm_multiplier_VI);
    }
    if (output_dset->dataspace_norm_multiplier_I >= 0) {
      H5Sclose(output_dset->dataspace_norm_multiplier_I);
    }
    if (output_dset->dataspace_norm_multiplier_QI >= 0) {
      H5Sclose(output_dset->dataspace_norm_multiplier_QI);
    }
    if (output_dset->dataspace_norm_multiplier_UI >= 0) {
      H5Sclose(output_dset->dataspace_norm_multiplier_UI);
    }
    if (output_dset->dataspace_norm_multiplier_VI >= 0) {
      H5Sclose(output_dset->dataspace_norm_multiplier_VI);
    }

    if (output_dset->datatype_id >= 0) {
      H5Tclose(output_dset->datatype_id);
    }
    if (output_dset->datatype_f32_id >= 0) {
      H5Tclose(output_dset->datatype_f32_id);
    }

    if (output_dset->group_id >= 0) {
      H5Gclose(output_dset->group_id);
    }
    output_dset->is_open = false;
  }

  free(output_dset);
}

/////////////////////////////////////////////////////
// THDF_write_3D_field_dataset_to_hdf5
/////////////////////////////////////////////////////
int
THDF_write_3D_field_dataset_to_hdf5(THDF_3D_field_handler_t *output_dset,    //
                                    THDF_3D_field_t         *output_field,   //
                                    hsize_t                  start_i,        //
                                    hsize_t                  start_j,        //
                                    hsize_t                  start_k,        //
                                    hsize_t                  start_incl,     //
                                    hsize_t                  start_azimuth,  //
                                    hsize_t                  count_i,        //
                                    hsize_t                  count_j,        //
                                    hsize_t                  count_k,        //
                                    hsize_t                  count_incl,     //
                                    hsize_t                  count_azimuth,  //
                                    hsize_t                  count_frequencies) {             //

  if (!output_dset || !output_field || !output_dset->is_open) {
    fprintf(stderr, "Invalid output dataset handler or field data for writing 3D field dataset\n");
    return -1;
  }

  hsize_t start[6] = {start_i, start_j, start_k, start_incl, start_azimuth, 0};
  hsize_t count[6] = {count_i, count_j, count_k, count_incl, count_azimuth, count_frequencies};

  // Memory space should match data layout: [x, y, z, incl, azimuth, frequencies]
  hid_t memspace = H5Screate_simple(6, count, NULL);

  // Write Stokes I, Q, U, V
  write_stokes_hyperslab(output_dset->dataset_id_I, output_dset->datatype_f32_id, memspace, start, count,
                         output_field->stokes_I, "Stokes I");
  write_stokes_hyperslab(output_dset->dataset_id_Q, output_dset->datatype_f32_id, memspace, start, count,
                         output_field->stokes_QI, "Stokes QI");
  write_stokes_hyperslab(output_dset->dataset_id_U, output_dset->datatype_f32_id, memspace, start, count,
                         output_field->stokes_UI, "Stokes UI");
  write_stokes_hyperslab(output_dset->dataset_id_V, output_dset->datatype_f32_id, memspace, start, count,
                         output_field->stokes_VI, "Stokes VI");

  H5Sclose(memspace);

  if (output_dset->normalized_output) {
    // Write normalization multipliers if normalized output is enabled
    hsize_t norm_start[5] = {start_i, start_j, start_k, start_incl, start_azimuth};
    hsize_t norm_count[5] = {count_i, count_j, count_k, count_incl, count_azimuth};
    hid_t   norm_memspace = H5Screate_simple(5, norm_count, NULL);

    write_stokes_hyperslab(output_dset->dataset_norm_multiplier_I, output_dset->datatype_id, norm_memspace, norm_start,
                           norm_count, output_field->norm_multiplier_I, "Normalization multiplier I");
    write_stokes_hyperslab(output_dset->dataset_norm_multiplier_QI, output_dset->datatype_id, norm_memspace, norm_start,
                           norm_count, output_field->norm_multiplier_QI, "Normalization multiplier QI");
    write_stokes_hyperslab(output_dset->dataset_norm_multiplier_UI, output_dset->datatype_id, norm_memspace, norm_start,
                           norm_count, output_field->norm_multiplier_UI, "Normalization multiplier UI");
    write_stokes_hyperslab(output_dset->dataset_norm_multiplier_VI, output_dset->datatype_id, norm_memspace, norm_start,
                           norm_count, output_field->norm_multiplier_VI, "Normalization multiplier VI");

    H5Sclose(norm_memspace);
  }

  return 0;
}

/////////////////////////////////////////////////////
// THDF_copy_3D_block_field
/////////////////////////////////////////////////////
void
THDF_copy_3D_block_field(THDF_3D_field_t *field,              //
                         const bool       normalized_output,  //
                         THDF_float_t    *stokes_IQUI,        // Input data in float 64 bits
                         THDF_float32_t  *stokes_out_I,       // Output normalized buffer data in float 32 bits
                         THDF_float32_t  *stokes_out_QI,      // This is a preallocated buffer that
                         THDF_float32_t  *stokes_out_UI,      // will hold the normalized Stokes Q/I values
                         THDF_float32_t  *stokes_out_VI,      //
                         THDF_float_t *norm_multiplier_I,  // Preallocated buffers that will hold the output normalization
                         THDF_float_t *norm_multiplier_QI,  // multipliers for each Stokes parameter
                         THDF_float_t *norm_multiplier_UI,  // Set to NULL if normalized_output is false
                         THDF_float_t *norm_multiplier_VI,  //
                         const int     N0_in_local_x,       //
                         const int     N0_in_local_y,       //
                         const int     N0_in_local_z,       //
                         const int     N_local_x,           //
                         const int     N_local_y,           //
                         const int     N_local_z,           //
                         const int     N_incl,              //
                         const int     N_azimuth,           //
                         const int     N_frequencies,       //
                         const int stride_in_x,  // Strides for the input data layout, allowing for flexible arrangements
                         const int stride_in_y,  // The function will compute the correct input index using
                         const int stride_in_z,  // these strides to access the data in the input arrays
                         const int stride_in_incl,         //
                         const int stride_in_azimuth,      //
                         const int stride_in_frequencies,  //
                         const int stride_in_stokes) {     //

  // Attention the data allocation is fully managed by the caller, this function only fills the provided buffers with
  // normalized data and normalization multipliers.

  //   const int total_local_points = N_local_x * N_local_y * N_local_z * N_incl * N_azimuth * N_frequencies;
  //   const int total_norm_points  = N_local_x * N_local_y * N_local_z * N_incl * N_azimuth;

  if (normalized_output) {
    field->norm_multiplier_I  = norm_multiplier_I;
    field->norm_multiplier_QI = norm_multiplier_QI;
    field->norm_multiplier_UI = norm_multiplier_UI;
    field->norm_multiplier_VI = norm_multiplier_VI;
  } else {
    field->norm_multiplier_I  = NULL;
    field->norm_multiplier_QI = NULL;
    field->norm_multiplier_UI = NULL;
    field->norm_multiplier_VI = NULL;
  }

  field->stokes_I  = stokes_out_I;
  field->stokes_QI = stokes_out_QI;
  field->stokes_UI = stokes_out_UI;
  field->stokes_VI = stokes_out_VI;

  for (int i = 0; i < N_local_x; i++) {
    for (int j = 0; j < N_local_y; j++) {
      for (int k = 0; k < N_local_z; k++) {
        for (int incl = 0; incl < N_incl; incl++) {
          for (int az = 0; az < N_azimuth; az++) {

            int norm_idx = ((((i * N_local_y + j) * N_local_z + k) * N_incl + incl) * N_azimuth) + az;

            double abs_max_I  = 0.0;
            double abs_max_QI = 0.0;
            double abs_max_UI = 0.0;
            double abs_max_VI = 0.0;

            if (normalized_output) {
              for (int f = 0; f < N_frequencies; f++) {
                const int in_index_I =
                    get_input_index(i, j, k, incl, az, f, 0, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
                                    stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
                                    stride_in_stokes);  // Stokes I is at offset 0 in the stokes dimension

                const int in_index_QI =
                    get_input_index(i, j, k, incl, az, f, 1, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
                                    stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
                                    stride_in_stokes);  // Stokes Q/I is at offset 1 in the stokes dimension

                const int in_index_UI =
                    get_input_index(i, j, k, incl, az, f, 2, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
                                    stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
                                    stride_in_stokes);  // Stokes U/I is at offset 2 in the stokes dimension

                const int in_index_VI =
                    get_input_index(i, j, k, incl, az, f, 3, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
                                    stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
                                    stride_in_stokes);  // Stokes V/I is at offset 3 in the stokes dimension

                if (fabs(stokes_IQUI[in_index_I]) > abs_max_I) {
                  abs_max_I = fabs(stokes_IQUI[in_index_I]);
                }

                if (fabs(stokes_IQUI[in_index_QI]) > abs_max_QI) {
                  abs_max_QI = fabs(stokes_IQUI[in_index_QI]);
                }

                if (fabs(stokes_IQUI[in_index_UI]) > abs_max_UI) {
                  abs_max_UI = fabs(stokes_IQUI[in_index_UI]);
                }

                if (fabs(stokes_IQUI[in_index_VI]) > abs_max_VI) {
                  abs_max_VI = fabs(stokes_IQUI[in_index_VI]);
                }
              }
            }

            const double norm_multiplier_I_  = (abs_max_I > 0.0) ? abs_max_I : 1.0;
            const double norm_multiplier_QI_ = (abs_max_QI > 0.0) ? abs_max_QI : 1.0;
            const double norm_multiplier_UI_ = (abs_max_UI > 0.0) ? abs_max_UI : 1.0;
            const double norm_multiplier_VI_ = (abs_max_VI > 0.0) ? abs_max_VI : 1.0;

            const double inv_norm_multiplier_I_  = 1.0 / norm_multiplier_I_;
            const double inv_norm_multiplier_QI_ = 1.0 / norm_multiplier_QI_;
            const double inv_norm_multiplier_UI_ = 1.0 / norm_multiplier_UI_;
            const double inv_norm_multiplier_VI_ = 1.0 / norm_multiplier_VI_;

            if (normalized_output) {
              norm_multiplier_I[norm_idx]  = norm_multiplier_I_;
              norm_multiplier_QI[norm_idx] = norm_multiplier_QI_;
              norm_multiplier_UI[norm_idx] = norm_multiplier_UI_;
              norm_multiplier_VI[norm_idx] = norm_multiplier_VI_;
            }

            for (int f = 0; f < N_frequencies; f++) {
              int idx = ((((((i)*N_local_y + (j)) * N_local_z + (k)) *  //
                               N_incl +
                           incl) *
                              N_azimuth +
                          az) *
                         N_frequencies) +
                        f;

              const int in_index_I =
                  get_input_index(i, j, k, incl, az, f, 0, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
                                  stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
                                  stride_in_stokes);  // Stokes I is at offset 0 in the stokes dimension

              const int in_index_QI =
                  get_input_index(i, j, k, incl, az, f, 1, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
                                  stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
                                  stride_in_stokes);  // Stokes Q/I is at offset 1 in the stokes dimension

              const int in_index_UI =
                  get_input_index(i, j, k, incl, az, f, 2, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
                                  stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
                                  stride_in_stokes);  // Stokes U/I is at offset 2 in the stokes dimension

              const int in_index_VI =
                  get_input_index(i, j, k, incl, az, f, 3, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
                                  stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
                                  stride_in_stokes);  // Stokes V/I is at offset 3 in the stokes dimension

              field->stokes_I[idx] = (THDF_float32_t)((THDF_float_t)stokes_IQUI[in_index_I] *  //
                                                      (THDF_float_t)inv_norm_multiplier_I_);   //

              field->stokes_QI[idx] = (THDF_float32_t)((THDF_float_t)stokes_IQUI[in_index_QI] *  //
                                                       (THDF_float_t)inv_norm_multiplier_QI_);   //

              field->stokes_UI[idx] = (THDF_float32_t)((THDF_float_t)stokes_IQUI[in_index_UI] *  //
                                                       (THDF_float_t)inv_norm_multiplier_UI_);   //

              field->stokes_VI[idx] = (THDF_float32_t)((THDF_float_t)stokes_IQUI[in_index_VI] *  //
                                                       (THDF_float_t)inv_norm_multiplier_VI_);   //
            }
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////
// write_3d_field_block_to_hdf5
//////////////////////////////////////////////////////
static void
write_3d_field_block_to_hdf5(const char *filename, MPI_Comm comm, bool normalized_output, int N_x, int N_y, int N_z,
                             int N_incl, int N_azimuth, int N_frequencies, int N_local_x, int N_local_y, int N_local_z,
                             int local_start_x, int local_start_y, int local_start_z, THDF_3D_field_t *output_field) {

  hid_t file_id = THDF_open_file_MPI(filename, comm);

  if (file_id < 0) {
    MPI_Abort(comm, -1);
  }

  int comm_rank;
  MPI_Comm_rank(comm, &comm_rank);
  const double write_start = MPI_Wtime();

  THDF_3D_field_handler_t *field_handler = THDF_create_3D_field_handler_mpi(file_id,                            //
                                                                            normalized_output,                  //
                                                                            N_x, N_y, N_z,                      //
                                                                            N_incl, N_azimuth, N_frequencies);  //

  THDF_write_3D_field_dataset_to_hdf5(field_handler, output_field, local_start_x, local_start_y, local_start_z, 0, 0,
                                      N_local_x, N_local_y, N_local_z, N_incl, N_azimuth, N_frequencies);

  const double write_elapsed = MPI_Wtime() - write_start;
  if (comm_rank == 0) {
    printf("MPI write elapsed time: %.6f s\n", write_elapsed);
  }

  THDF_close_3D_field_handler_mpi(field_handler);
  H5Fclose(file_id);
}

//////////////////////////////////////////////////////
// write_3d_field_block_mpi
//////////////////////////////////////////////////////
void
write_3d_field_block_mpi(const char     *filename,            //
                         MPI_Comm        comm,                //
                         bool            normalized_output,   //
                         int             N_x,                 //
                         int             N_y,                 //
                         int             N_z,                 //
                         int             N_incl,              //
                         int             N_azimuth,           //
                         int             N_frequencies,       //
                         int             N_stokes,            //
                         int             N_local_x,           //
                         int             N_local_y,           //
                         int             N_local_z,           //
                         int             local_start_x,       //
                         int             local_start_y,       //
                         int             local_start_z,       //
                         double         *stokes_IQUI,         //
                         THDF_float32_t *stokes_I,            //
                         THDF_float32_t *stokes_QI,           //
                         THDF_float32_t *stokes_UI,           //
                         THDF_float32_t *stokes_VI,           //
                         THDF_float_t   *norm_multiplier_I,   //
                         THDF_float_t   *norm_multiplier_QI,  //
                         THDF_float_t   *norm_multiplier_UI,  //
                         THDF_float_t   *norm_multiplier_VI,  //
                         int             stride_x,            //
                         int             stride_y,            //
                         int             stride_z,            //
                         int             stride_incl,         //
                         int             stride_azimuth,      //
                         int             stride_frequencies,  //
                         int             stride_stokes) {                 //
                                                              //
  (void)N_stokes;

  THDF_3D_field_t output_field;

  THDF_copy_3D_block_field(&output_field, normalized_output, stokes_IQUI, stokes_I, stokes_QI, stokes_UI, stokes_VI,
                           norm_multiplier_I, norm_multiplier_QI, norm_multiplier_UI, norm_multiplier_VI, 0, 0, 0,
                           N_local_x, N_local_y, N_local_z, N_incl, N_azimuth, N_frequencies, stride_x, stride_y,
                           stride_z, stride_incl, stride_azimuth, stride_frequencies, stride_stokes);

  write_3d_field_block_to_hdf5(filename, comm, normalized_output, N_x, N_y, N_z, N_incl, N_azimuth, N_frequencies,
                               N_local_x, N_local_y, N_local_z, local_start_x, local_start_y, local_start_z,
                               &output_field);
}