#ifndef __THDF_OUTPUT_3D_FIELD_H__
#define __THDF_OUTPUT_3D_FIELD_H__

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf_output_field.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {

  // The data must be stored in a ROW MAJOR 5D contiguous array:
  // with memory layout match data layout: [x, y, incl, azimuth, frequencies]
  THDF_float32_t *stokes_I;   //
  THDF_float32_t *stokes_QI;  //
  THDF_float32_t *stokes_UI;  //
  THDF_float32_t *stokes_VI;  //

  THDF_n_float_t *norm_multiplier_I;   //
  THDF_n_float_t *norm_multiplier_QI;  // Normalization multipliers for each Stokes parameter
  THDF_n_float_t *norm_multiplier_UI;  // At the moment not used
  THDF_n_float_t *norm_multiplier_VI;  //

} THDF_3D_field_t;

typedef struct {
  hid_t file_id;
  hid_t group_id;

  hid_t dataset_id_I;
  hid_t dataset_id_Q;
  hid_t dataset_id_U;
  hid_t dataset_id_V;

  hid_t dataspace_id_I;
  hid_t dataspace_id_Q;
  hid_t dataspace_id_U;
  hid_t dataspace_id_V;

  hid_t dataset_norm_multiplier_I;
  hid_t dataset_norm_multiplier_QI;
  hid_t dataset_norm_multiplier_UI;
  hid_t dataset_norm_multiplier_VI;

  hid_t dataspace_norm_multiplier_I;
  hid_t dataspace_norm_multiplier_QI;
  hid_t dataspace_norm_multiplier_UI;
  hid_t dataspace_norm_multiplier_VI;

  hid_t datatype_id;
  hid_t datatype_f32_id;

  int N_x;
  int N_y;
  int N_z;
  int N_incl;
  int N_azimuth;
  int N_frequencies;

  bool is_open;
  bool normalized_output;  // Flag to indicate if the output is normalized or not
} THDF_3D_field_handler_t;

int
THDF_create_3D_field_handler_mpi(const char *path, const bool normalized_output, const int N_x, const int N_y,
                                 const int N_z, const int N_incl, const int N_azimuth, const int N_frequencies);

THDF_3D_field_handler_t *
THDF_open_3D_field_handler_mpi(hid_t file, const bool normalized_output, const int N_x, const int N_y, const int N_z,
                               const int N_incl, const int N_azimuth, const int N_frequencies);

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
                                    hsize_t                  count_frequencies);

/**
 * @brief Close and free resources associated with an HDF5 3D field handler (MPI version).
 *
 * @param output_dset
 */
void
THDF_close_3D_field_handler_mpi(THDF_3D_field_handler_t *output_dset);

/**
 * @brief Copy and reshuffle a 3D block field using local dimensions.
 * The input data is expected to be in a 1D array with the following memory layout: [x, y, z, incl, azimuth, frequencies]
 * The output data is reshuffled to match the HDF5 dataset layout: [x, y, z, incl, azimuth, frequencies] but stored in a
 * 1D array in row-major order. The normalization is applied if the normalized_output flag is set to true, in this case
 * the function also computes the normalization multipliers for each Stokes parameter and fills the provided buffers.
 * The stride parameters indicate the memory layout of the input data, allowing for flexible data arrangements.
 *
 * @param field
 * @param N_local_x
 * @param N_local_y
 * @param N_local_z
 * @param N_incl
 * @param N_azimuth
 * @param N_frequencies
 */
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
                         const int stride_in_stokes);      //

void
write_3d_field_block_mpi(hid_t           file_id,             //
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
                         int             stride_stokes);      //

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // __THDF_OUTPUT_3D_FIELD_H__