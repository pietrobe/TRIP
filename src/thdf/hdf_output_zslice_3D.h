#ifndef HDF_OUTPUT_ZSLICE_3D_H
#define HDF_OUTPUT_ZSLICE_3D_H

#if defined(__cplusplus)
extern "C" {
#endif

#include <hdf5.h>
#include <mpi.h>
#include <stdbool.h>

#include "hdf_output_3D_field.h"

typedef struct {
  hid_t file_id;
  hid_t group_id;
  hid_t dataset_id_I, dataset_id_Q, dataset_id_U, dataset_id_V;
  hid_t dataset_norm_multiplier_I, dataset_norm_multiplier_QI;
  hid_t dataset_norm_multiplier_UI, dataset_norm_multiplier_VI;
  hid_t datatype_f32_id;
  hid_t datatype_f64_id;
  bool  is_open;
  bool  normalized_output;
  int   N_x, N_y, N_z, N_incl, N_azimuth, N_frequencies;
} THDF_zslab_handler_t;

typedef struct {
  MPI_Comm slab_comm;
  MPI_Comm root_comm;
  int      slab_id;
  int      slab_rank;
  int      slab_size;
  int      n_slabs;
  int      N_local_z;
} THDF_zslab_topology_t;

THDF_zslab_topology_t
build_zslab_topology(MPI_Comm world, int local_start_z, int N_local_z, int mpi_rank);

void
build_zslab_filename(char *buf, size_t sz, const char *dir, int slab_id);

void
build_slab_filename_slice(char *buf, size_t sz,                        //
                          const int glob_slice_start_z, int slice_nz,  //
                          const int glob_slice_start_x, int slice_nx,  //
                          const int glob_slice_start_y, int slice_ny);

int
THDF_create_zslab_file(MPI_Comm comm, const char *path, bool normalized_output, int N_x, int N_y, int N_local_z,
                       int N_incl, int N_azimuth, int N_frequencies);

THDF_zslab_handler_t *
THDF_open_zslab_handler(hid_t file, bool normalized_output, int N_x, int N_y, int N_local_z, int N_incl, int N_azimuth,
                        int N_frequencies);

void
THDF_close_zslab_handler(THDF_zslab_handler_t *h);

int
THDF_write_zslab_block(THDF_zslab_handler_t *h, const THDF_3D_field_t *field, int local_start_x, int local_start_y,
                       int local_zi, int N_local_x, int N_local_y, int current_nz, int N_incl, int N_azimuth,
                       int N_frequencies);

void
create_zslab_vds_master(const char *dir, const char *master_name, int n_slabs, int N_x, int N_y, int N_z, int N_incl,
                        int N_azimuth, int N_freq, bool normalized_output, const int *slab_z_start, const int *slab_nz);

#if defined(__cplusplus)
}
#endif

#endif