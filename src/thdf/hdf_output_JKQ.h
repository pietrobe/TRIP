#ifndef __HDF_OUTPUT_JKQ_H__
#define __HDF_OUTPUT_JKQ_H__

#include "hdf_output_field.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef float  THDF_JKQ_float_t;
typedef double THDF_JKQ_double_t;
typedef double THDF_JKQ_n_float_t;

/**
 * @brief Structure to hold the Stokes parameters field data for arbitrary directions
 * at a given spatial location that will be written to or read from HDF5 files.
 */
typedef struct {

  THDF_JKQ_float_t *JKQ_re;  // mean intensity real part for PRD
  THDF_JKQ_float_t *JKQ_im;  // mean intensity imaginary part for PRD

  THDF_JKQ_n_float_t *norm_multiplier_JKQ_re;  // Normalization multipliers for each Stokes parameter
  THDF_JKQ_n_float_t *norm_multiplier_JKQ_im;  // Normalization multipliers for each Stokes parameter

} THDF_JKQ_field_t;

/**
 * @brief Structure to handle HDF5 output field dataset operations.
 */
typedef struct {
  hid_t file_id;
  hid_t group_id;

  hid_t dataset_id_JKQ_re;
  hid_t dataset_id_JKQ_im;

  hid_t dataspace_id_JKQ_re;
  hid_t dataspace_id_JKQ_im;

  hid_t dataset_JKQ_norm_mult_re;
  hid_t dataset_JKQ_norm_mult_im;

  hid_t dataspace_JKQ_norm_mult_re;
  hid_t dataspace_JKQ_norm_mult_im;

  hid_t datatype_id;
  hid_t datatype_norm_mult_id;

  int N_x;
  int N_y;
  int N_z;
  int N_JKQ_real;
  int N_JKQ_imag;
  int N_frequencies;

  bool is_open;
} THDF_JKQ_handler_t;

/**
 * @brief Structure to hold the KQ matrix used in the JKQ formalism.
 */
typedef struct {

  int KQ_size;  // The number of elements (rows) in the KQ matrix
                // note: the number of columns is always 2 (K and Q)

  int *KQ_table;  // The KQ matrix data stored in a 1D array (row-major order)

  ////////////////////////////////////////////////////////////////////////////////////
  int KQ_compressed_size_real;  // Size of compressed indices for real part
  int KQ_compressed_size_imag;  // Size of compressed indices for imaginary part

  int *KQ_compressed_real;  // Compressed indices for real part
  int *KQ_compressed_imag;  // Compressed indices for imaginary part

} THDF_KQ_table_t;

/**
 * @brief Create and initialize a JKQ field handler for HDF5 output with MPI support.
 * @param file HDF5 file identifier.
 * @param N_x Number of grid points in the x-dimension.
 * @param N_y Number of grid points in the y-dimension.
 * @param N_z Number of grid points in the z-dimension.
 * @param N_JKQ_real Number of JKQ components for the real part.
 * @param N_JKQ_imag Number of JKQ components for the imaginary part.
 * @param N_frequencies Number of frequency points.
 * @return Pointer to the initialized JKQ field handler, or NULL on failure.
 */
THDF_JKQ_handler_t *
THDF_create_JKQ_field_handler_mpi(hid_t file,            //
                                  int   N_x,             //
                                  int   N_y,             //
                                  int   N_z,             //
                                  int   N_JKQ_real,      //
                                  int   N_JKQ_imag,      //
                                  int   N_frequencies);  //

/**
 * @brief Close and free resources associated with the JKQ field handler.
 * @param handler Pointer to the JKQ field handler to close.
 */
void
THDF_close_JKQ_field_handler_mpi(THDF_JKQ_handler_t *handler);  //

/**
 * @brief Write the JKQ field data to the HDF5 file.
 * @param handler Pointer to the JKQ field handler.
 * @param JKQ_field Pointer to the JKQ field data to write.
 * @param start_i Starting index in the x-dimension.
 * @param start_j Starting index in the y-dimension.
 * @param start_k Starting index in the z-dimension.
 * @param count_i Number of elements to write in the x-dimension.
 * @param count_j Number of elements to write in the y-dimension.
 * @param count_k Number of elements to write in the z-dimension.
 * @return 0 on success, non-zero on failure.
 */
int
THDF_write_JKQ_field_to_hdf5(THDF_JKQ_handler_t *handler,    //
                             THDF_JKQ_field_t   *JKQ_field,  //
                             hsize_t             start_i,    //
                             hsize_t             start_j,    //
                             hsize_t             start_k,    //
                             hsize_t             count_i,    //
                             hsize_t             count_j,    //
                             hsize_t             count_k);   //

int
THDF_write_KQ_table_to_hdf5(hid_t            file,       //
                            THDF_KQ_table_t *KQ_table);  //

#ifdef __cplusplus
}
#endif

#endif  // __§HDF_OUTPUT_JKQ_H__