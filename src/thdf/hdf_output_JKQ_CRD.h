#ifndef __HDF_OUTPUT_JKQ_CRD_H__
#define __HDF_OUTPUT_JKQ_CRD_H__

#include "hdf_output_JKQ.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {

  THDF_JKQ_double_t *JKQ_re;  // Mean intensity real part for CRD.
  THDF_JKQ_double_t *JKQ_im;  // Mean intensity imaginary part for CRD.

} THDF_JKQ_CRD_field_t;

/*
 * Create (or open) HDF5 datasets /JKQ/JKQ_re and /JKQ/JKQ_im for CRD output.
 *
 * The global dataset layout is:
 *   JKQ_re[N_x][N_y][N_z][N_JKQ_real]
 *   JKQ_im[N_x][N_y][N_z][N_JKQ_imag]
 *
 * Returns a valid handler on success, NULL on failure.
 */

THDF_JKQ_handler_t *
THDF_create_JKQ_CRD_field_handler_mpi(hid_t file,         //
                                      int   N_x,          //
                                      int   N_y,          //
                                      int   N_z,          //
                                      int   N_JKQ_real,   //
                                      int   N_JKQ_imag);  //

/*
 * Write a local 3D subdomain of JKQ CRD values using hyperslabs.
 *
 * The written block starts at (start_i, start_j, start_k) in the global
 * [N_x, N_y, N_z] domain and spans (count_i, count_j, count_k).
 *
 * Buffer expectations:
 *   JKQ_field->JKQ_re has count_i*count_j*count_k*N_JKQ_real elements.
 *   JKQ_field->JKQ_im has count_i*count_j*count_k*N_JKQ_imag elements.
 *
 * Returns 0 on success, -1 on failure.
 */

int
THDF_write_JKQ_CRD_field_to_hdf5(THDF_JKQ_handler_t   *handler,    //
                                 THDF_JKQ_CRD_field_t *JKQ_field,  //
                                 hsize_t               start_i,    //
                                 hsize_t               start_j,    //
                                 hsize_t               start_k,    //
                                 hsize_t               count_i,    //
                                 hsize_t               count_j,    //
                                 hsize_t               count_k);   //

/*
 * Close HDF5 resources owned by the CRD handler and free the handler.
 * Safe to call with NULL.
 */

void
THDF_close_JKQ_CRD_field_handler_mpi(THDF_JKQ_handler_t *handler);  //

/*
 * Write KQ metadata for CRD output.
 * This forwards to the shared KQ table writer.
 *
 * Returns 0 on success, negative value on failure.
 */

int
THDF_write_KQ_CRD_table_to_hdf5(hid_t            file,       //
                                THDF_KQ_table_t *KQ_table);  //

#ifdef __cplusplus
}
#endif

#endif  // __HDF_OUTPUT_JKQ_CRD_H__