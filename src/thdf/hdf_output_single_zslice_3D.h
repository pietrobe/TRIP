#ifndef HDF_OUTPUT_SINGLE_ZSLICE_3D_H
#define HDF_OUTPUT_SINGLE_ZSLICE_3D_H

#include <hdf5.h>
#include <mpi.h>
#include <stdbool.h>

#include "hdf_output_zslice_3D.h"

#ifdef __cplusplus
extern "C" {
#endif

int
THDF_create_zslab_file_single(MPI_Comm comm, const char *path, bool normalized_output, int N_x, int N_y, int N_z,
                              int N_incl, int N_azimuth, int N_frequencies);

THDF_zslab_handler_t *
THDF_open_zslab_handler_single(hid_t file, bool normalized_output, int N_x, int N_y, int N_z_global, int N_incl,
                               int N_azimuth, int N_frequencies);

int
THDF_write_zslab_block_single(THDF_zslab_handler_t *handler, const THDF_3D_field_t *field, int local_start_x,
                              int local_start_y, int global_start_z, int N_local_x, int N_local_y, int current_nz,
                              int N_incl, int N_azimuth, int N_frequencies);

void
THDF_close_zslab_handler_single(THDF_zslab_handler_t *h);

#ifdef __cplusplus
}
#endif

#endif