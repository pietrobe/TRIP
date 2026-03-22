#include "hdf_output_zslice_3D.h"
#include "hdf_fapl_from_env.h"

#include <hdf5.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#define HDF_CLOSE_IF(id, fn) \
  do {                       \
    if ((id) >= 0) {         \
      (fn)(id);              \
      (id) = -1;             \
    }                        \
  } while (0)

#define THDF_CHECK(call, msg, label)      \
  do {                                    \
    if ((call) < 0) {                     \
      fprintf(stderr, "[HDF] " msg "\n"); \
      rc = -1;                            \
      goto label;                         \
    }                                     \
  } while (0)

THDF_zslab_topology_t
build_zslab_topology(MPI_Comm world, int local_start_z, int N_local_z, int mpi_rank) {
  THDF_zslab_topology_t t;
  t.slab_id   = local_start_z;
  t.N_local_z = N_local_z;
  t.root_comm = MPI_COMM_NULL;

  MPI_Comm_split(world, local_start_z, mpi_rank, &t.slab_comm);
  MPI_Comm_rank(t.slab_comm, &t.slab_rank);
  MPI_Comm_size(t.slab_comm, &t.slab_size);

  int is_root = (t.slab_rank == 0) ? 0 : MPI_UNDEFINED;
  MPI_Comm_split(world, is_root, mpi_rank, &t.root_comm);

  if (t.slab_rank == 0)
    MPI_Comm_size(t.root_comm, &t.n_slabs);

  MPI_Bcast(&t.n_slabs, 1, MPI_INT, 0, t.slab_comm);

  return t;
}

void
build_zslab_filename(char *buf, size_t sz, const char *dir, int slab_id) {
  snprintf(buf, sz, "%s/slab_z%05d.h5", dir, slab_id);
}

void
build_slab_filename_slice(char *buf, size_t sz,                        //
                          const int glob_slice_start_z, int slice_nz,  //
                          const int glob_slice_start_x, int slice_nx,  //
                          const int glob_slice_start_y, int slice_ny)  //
{

  const int end_slice_z = glob_slice_start_z + slice_nz;
  const int end_slice_x = glob_slice_start_x + slice_nx;
  const int end_slice_y = glob_slice_start_y + slice_ny;

  snprintf(buf, sz, "slice_x%d-%d_y%d-%d_z%d-%d.h5", glob_slice_start_x, end_slice_x, glob_slice_start_y, end_slice_y,
           glob_slice_start_z, end_slice_z);
}

int
THDF_create_zslab_file(MPI_Comm comm, const char *path, bool normalized_output, int N_x, int N_y, int N_local_z,
                       int N_incl, int N_azimuth, int N_frequencies) {
  hid_t fapl = -1, file = -1, group = -1;
  hid_t dcpl = -1, dcpl_n = -1;
  hid_t fspace = -1, nspace = -1;
  hid_t f32 = -1, f64 = -1;
  hid_t dset = -1;
  int   rc   = 0;
  int   rank;
  MPI_Comm_rank(comm, &rank);

  fapl = build_fapl_from_env(comm);
  if (fapl < 0)
    return -1;

  file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  HDF_CLOSE_IF(fapl, H5Pclose); /* chiudi subito dopo H5Fcreate */

  group = H5Gcreate2(file, "/radiation_field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (group < 0) {
    rc = -1;
    goto cleanup;
  }

  f32 = H5Tcopy(THDF_get_hdf_float32_datatype());
  if (f32 < 0) {
    rc = -1;
    goto cleanup;
  }

  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  if (dcpl < 0) {
    rc = -1;
    goto cleanup;
  }
  {
    hsize_t cy        = (hsize_t)(N_y < 16 ? N_y : 16);
    hsize_t cz        = (hsize_t)(N_local_z < 16 ? N_local_z : 16);
    hsize_t chunk6[6] = {1, cy, cz, (hsize_t)N_incl, (hsize_t)N_azimuth, (hsize_t)N_frequencies};
    THDF_CHECK(H5Pset_chunk(dcpl, 6, chunk6), "chunk 6D", cleanup);
  }

  {
    hsize_t dims6[6] = {(hsize_t)N_x,    (hsize_t)N_y,       (hsize_t)N_local_z,
                        (hsize_t)N_incl, (hsize_t)N_azimuth, (hsize_t)N_frequencies};
    fspace           = H5Screate_simple(6, dims6, NULL);
    if (fspace < 0) {
      rc = -1;
      goto cleanup;
    }

    const char *names6[4] = {"I", "QI_pc", "UI_pc", "VI_pc"};
    for (int i = 0; i < 4; i++) {
      dset = H5Dcreate2(group, names6[i], f32, fspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);
      if (dset < 0) {
        if (rank == 0)
          fprintf(stderr, "[create] Error creating dataset %s\n", names6[i]);
        rc = -1;
        goto cleanup;
      }
      HDF_CLOSE_IF(dset, H5Dclose);
    }
    HDF_CLOSE_IF(fspace, H5Sclose);
  }
  HDF_CLOSE_IF(dcpl, H5Pclose);

  if (!normalized_output)
    goto sync;

  f64 = H5Tcopy(THDF_get_hdf_float_datatype());
  if (f64 < 0) {
    rc = -1;
    goto cleanup;
  }

  dcpl_n = H5Pcreate(H5P_DATASET_CREATE);
  if (dcpl_n < 0) {
    rc = -1;
    goto cleanup;
  }
  {
    hsize_t cy        = (hsize_t)(N_y < 16 ? N_y : 16);
    hsize_t cz        = (hsize_t)(N_local_z < 16 ? N_local_z : 16);
    hsize_t chunk5[5] = {1, cy, cz, (hsize_t)N_incl, (hsize_t)N_azimuth};
    THDF_CHECK(H5Pset_chunk(dcpl_n, 5, chunk5), "chunk 5D", cleanup);
  }

  {
    hsize_t dims5[5] = {(hsize_t)N_x, (hsize_t)N_y, (hsize_t)N_local_z, (hsize_t)N_incl, (hsize_t)N_azimuth};
    nspace           = H5Screate_simple(5, dims5, NULL);
    if (nspace < 0) {
      rc = -1;
      goto cleanup;
    }

    const char *names5[4] = {"norm_multiplier_I", "norm_multiplier_QI_pc", "norm_multiplier_UI_pc",
                             "norm_multiplier_VI_pc"};
    for (int i = 0; i < 4; i++) {
      dset = H5Dcreate2(group, names5[i], f64, nspace, H5P_DEFAULT, dcpl_n, H5P_DEFAULT);
      if (dset < 0) {
        if (rank == 0)
          fprintf(stderr, "[create] Error creating dataset %s\n", names5[i]);
        rc = -1;
        goto cleanup;
      }
      HDF_CLOSE_IF(dset, H5Dclose);
    }
    HDF_CLOSE_IF(nspace, H5Sclose);
  }
  HDF_CLOSE_IF(dcpl_n, H5Pclose);

sync:
  if (H5Fflush(file, H5F_SCOPE_GLOBAL) < 0) {
    if (rank == 0)
      fprintf(stderr, "[create] H5Fflush failed\n");
    rc = -1;
  }

cleanup:
  HDF_CLOSE_IF(dset, H5Dclose);
  HDF_CLOSE_IF(fspace, H5Sclose);
  HDF_CLOSE_IF(nspace, H5Sclose);
  HDF_CLOSE_IF(dcpl, H5Pclose);
  HDF_CLOSE_IF(dcpl_n, H5Pclose);
  HDF_CLOSE_IF(f32, H5Tclose);
  HDF_CLOSE_IF(f64, H5Tclose);
  HDF_CLOSE_IF(group, H5Gclose);
  HDF_CLOSE_IF(fapl, H5Pclose);
  if (file >= 0)
    H5Fclose(file);

  return rc;
}

THDF_zslab_handler_t *
THDF_open_zslab_handler(hid_t file, bool normalized_output, int N_x, int N_y, int N_local_z, int N_incl, int N_azimuth,
                        int N_frequencies) {

  THDF_zslab_handler_t *h = calloc(1, sizeof(*h));
  if (!h)
    return NULL;

  h->file_id                    = file;
  h->group_id                   = -1;
  h->dataset_id_I               = -1;
  h->dataset_id_Q               = -1;
  h->dataset_id_U               = -1;
  h->dataset_id_V               = -1;
  h->dataset_norm_multiplier_I  = -1;
  h->dataset_norm_multiplier_QI = -1;
  h->dataset_norm_multiplier_UI = -1;
  h->dataset_norm_multiplier_VI = -1;
  h->datatype_f32_id            = -1;
  h->datatype_f64_id            = -1;
  h->is_open                    = false;
  h->normalized_output          = normalized_output;
  h->N_x                        = N_x;
  h->N_y                        = N_y;
  h->N_z                        = N_local_z;
  h->N_incl                     = N_incl;
  h->N_azimuth                  = N_azimuth;
  h->N_frequencies              = N_frequencies;

  h->group_id = H5Gopen2(file, "/radiation_field", H5P_DEFAULT);
  if (h->group_id < 0)
    goto fail;

  h->datatype_f32_id = H5Tcopy(THDF_get_hdf_float32_datatype());
  if (h->datatype_f32_id < 0)
    goto fail;

  const char *names6[4] = {"I", "QI_pc", "UI_pc", "VI_pc"};
  hid_t      *dsets6[4] = {&h->dataset_id_I, &h->dataset_id_Q, &h->dataset_id_U, &h->dataset_id_V};
  for (int i = 0; i < 4; i++) {
    *dsets6[i] = H5Dopen2(h->group_id, names6[i], H5P_DEFAULT);
    if (*dsets6[i] < 0)
      goto fail;
  }

  if (!normalized_output) {
    h->is_open = true;
    return h;
  }

  h->datatype_f64_id = H5Tcopy(THDF_get_hdf_float_datatype());
  if (h->datatype_f64_id < 0)
    goto fail;

  const char *names5[4] = {"norm_multiplier_I", "norm_multiplier_QI_pc", "norm_multiplier_UI_pc",
                           "norm_multiplier_VI_pc"};
  hid_t      *dsets5[4] = {&h->dataset_norm_multiplier_I, &h->dataset_norm_multiplier_QI, &h->dataset_norm_multiplier_UI,
                           &h->dataset_norm_multiplier_VI};
  for (int i = 0; i < 4; i++) {
    *dsets5[i] = H5Dopen2(h->group_id, names5[i], H5P_DEFAULT);
    if (*dsets5[i] < 0)
      goto fail;
  }

  h->is_open = true;
  return h;

fail:
  HDF_CLOSE_IF(h->dataset_id_I, H5Dclose);
  HDF_CLOSE_IF(h->dataset_id_Q, H5Dclose);
  HDF_CLOSE_IF(h->dataset_id_U, H5Dclose);
  HDF_CLOSE_IF(h->dataset_id_V, H5Dclose);
  HDF_CLOSE_IF(h->dataset_norm_multiplier_I, H5Dclose);
  HDF_CLOSE_IF(h->dataset_norm_multiplier_QI, H5Dclose);
  HDF_CLOSE_IF(h->dataset_norm_multiplier_UI, H5Dclose);
  HDF_CLOSE_IF(h->dataset_norm_multiplier_VI, H5Dclose);
  HDF_CLOSE_IF(h->datatype_f32_id, H5Tclose);
  HDF_CLOSE_IF(h->datatype_f64_id, H5Tclose);
  HDF_CLOSE_IF(h->group_id, H5Gclose);
  free(h);
  return NULL;
}

void
THDF_close_zslab_handler(THDF_zslab_handler_t *h) {
  if (!h || !h->is_open)
    return;
  HDF_CLOSE_IF(h->dataset_id_I, H5Dclose);
  HDF_CLOSE_IF(h->dataset_id_Q, H5Dclose);
  HDF_CLOSE_IF(h->dataset_id_U, H5Dclose);
  HDF_CLOSE_IF(h->dataset_id_V, H5Dclose);
  HDF_CLOSE_IF(h->dataset_norm_multiplier_I, H5Dclose);
  HDF_CLOSE_IF(h->dataset_norm_multiplier_QI, H5Dclose);
  HDF_CLOSE_IF(h->dataset_norm_multiplier_UI, H5Dclose);
  HDF_CLOSE_IF(h->dataset_norm_multiplier_VI, H5Dclose);
  HDF_CLOSE_IF(h->datatype_f32_id, H5Tclose);
  HDF_CLOSE_IF(h->datatype_f64_id, H5Tclose);
  HDF_CLOSE_IF(h->group_id, H5Gclose);
  h->is_open = false;
  free(h);
}

static int
write_stokes_hyperslab(hid_t dataset_id, hid_t datatype_id, hid_t memspace, hid_t dxpl, const hsize_t *start,
                       const hsize_t *count, const void *data, const char *name, bool empty_io) {
  hid_t filespace = H5Dget_space(dataset_id);
  if (filespace < 0) {
    fprintf(stderr, "[write] H5Dget_space failed for %s\n", name);
    return -1;
  }

  herr_t sel_ret;
  if (empty_io) {
    sel_ret = H5Sselect_none(filespace);
  } else {
    sel_ret = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);
  }
  if (sel_ret < 0) {
    fprintf(stderr, "[write] H5S selection failed for %s\n", name);
    H5Sclose(filespace);
    return -1;
  }

  static const char dummy = 0;
  const void       *ptr   = data ? data : (const void *)&dummy;
  herr_t            ret   = H5Dwrite(dataset_id, datatype_id, memspace, filespace, dxpl, ptr);
  H5Sclose(filespace);

  if (ret < 0) {
    fprintf(stderr, "[write] H5Dwrite failed for %s\n", name);
    return -1;
  }
  return 0;
}

int
THDF_write_zslab_block(THDF_zslab_handler_t *h, const THDF_3D_field_t *field, int local_start_x, int local_start_y,
                       int local_zi, int N_local_x, int N_local_y, int current_nz, int N_incl, int N_azimuth,
                       int N_frequencies) {

  if (!h || !h->is_open)
    return -1;

  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  if (dxpl < 0)
    return -1;
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

  const bool empty_io = (current_nz <= 0);

  hsize_t start[6] = {(hsize_t)local_start_x, (hsize_t)local_start_y, (hsize_t)local_zi, 0, 0, 0};
  hsize_t count[6] = {(hsize_t)N_local_x, (hsize_t)N_local_y, (hsize_t)(empty_io ? 0 : current_nz),
                      (hsize_t)N_incl,    (hsize_t)N_azimuth, (hsize_t)N_frequencies};

  hid_t memspace = empty_io ? H5Screate(H5S_NULL) : H5Screate_simple(6, count, NULL);
  if (memspace < 0) {
    H5Pclose(dxpl);
    return -1;
  }

  int rc = 0;
  rc |= write_stokes_hyperslab(h->dataset_id_I, h->datatype_f32_id, memspace, dxpl, start, count, field->stokes_I,
                               "Stokes I", empty_io);
  rc |= write_stokes_hyperslab(h->dataset_id_Q, h->datatype_f32_id, memspace, dxpl, start, count, field->stokes_QI,
                               "Stokes Q", empty_io);
  rc |= write_stokes_hyperslab(h->dataset_id_U, h->datatype_f32_id, memspace, dxpl, start, count, field->stokes_UI,
                               "Stokes U", empty_io);
  rc |= write_stokes_hyperslab(h->dataset_id_V, h->datatype_f32_id, memspace, dxpl, start, count, field->stokes_VI,
                               "Stokes V", empty_io);
  H5Sclose(memspace);

#ifndef NDEBUG
  {
    uint32_t lc, gc;
    H5Pget_mpio_no_collective_cause(dxpl, &lc, &gc);
    if (gc != 0)
      fprintf(stderr, "[write] WARNING: collective I/O fallback cause=0x%x\n", gc);
  }
#endif

  if (!h->normalized_output) {
    H5Pclose(dxpl);
    return rc;
  }

  hsize_t nstart[5] = {(hsize_t)local_start_x, (hsize_t)local_start_y, (hsize_t)local_zi, 0, 0};
  hsize_t ncount[5] = {(hsize_t)N_local_x, (hsize_t)N_local_y, (hsize_t)(empty_io ? 0 : current_nz), (hsize_t)N_incl,
                       (hsize_t)N_azimuth};
  hid_t   nmemspace = empty_io ? H5Screate(H5S_NULL) : H5Screate_simple(5, ncount, NULL);
  if (nmemspace < 0) {
    H5Pclose(dxpl);
    return -1;
  }

  rc |= write_stokes_hyperslab(h->dataset_norm_multiplier_I, h->datatype_f64_id, nmemspace, dxpl, nstart, ncount,
                               field->norm_multiplier_I, "norm I", empty_io);
  rc |= write_stokes_hyperslab(h->dataset_norm_multiplier_QI, h->datatype_f64_id, nmemspace, dxpl, nstart, ncount,
                               field->norm_multiplier_QI, "norm Q", empty_io);
  rc |= write_stokes_hyperslab(h->dataset_norm_multiplier_UI, h->datatype_f64_id, nmemspace, dxpl, nstart, ncount,
                               field->norm_multiplier_UI, "norm U", empty_io);
  rc |= write_stokes_hyperslab(h->dataset_norm_multiplier_VI, h->datatype_f64_id, nmemspace, dxpl, nstart, ncount,
                               field->norm_multiplier_VI, "norm V", empty_io);
  H5Sclose(nmemspace);
  H5Pclose(dxpl);
  return rc;
}

void
create_zslab_vds_master(const char *dir, const char *master_name, int n_slabs, int N_x, int N_y, int N_z, int N_incl,
                        int N_azimuth, int N_freq, bool normalized_output, const int *slab_z_start, const int *slab_nz) {

  char master_path[PATH_MAX];
  snprintf(master_path, sizeof(master_path), "%s/%s", dir, master_name);

  hid_t master = H5Fcreate(master_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t grp    = H5Gcreate2(master, "/radiation_field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hsize_t vdims6[6] = {(hsize_t)N_x, (hsize_t)N_y, (hsize_t)N_z, (hsize_t)N_incl, (hsize_t)N_azimuth, (hsize_t)N_freq};
  hsize_t vdims5[5] = {(hsize_t)N_x, (hsize_t)N_y, (hsize_t)N_z, (hsize_t)N_incl, (hsize_t)N_azimuth};

  hid_t f32 = H5Tcopy(THDF_get_hdf_float32_datatype());
  hid_t f64 = H5Tcopy(THDF_get_hdf_float_datatype());

  const char *names6[4] = {"I", "QI_pc", "UI_pc", "VI_pc"};
  const char *names5[4] = {"norm_multiplier_I", "norm_multiplier_QI_pc", "norm_multiplier_UI_pc",
                           "norm_multiplier_VI_pc"};

  for (int di = 0; di < 4; di++) {
    hid_t vspace = H5Screate_simple(6, vdims6, NULL);
    hid_t dcpl   = H5Pcreate(H5P_DATASET_CREATE);

    for (int s = 0; s < n_slabs; s++) {
      char src[PATH_MAX];
      build_zslab_filename(src, sizeof(src), dir, slab_z_start[s]);

      hsize_t src_dims[6] = {(hsize_t)N_x,    (hsize_t)N_y,       (hsize_t)slab_nz[s],
                             (hsize_t)N_incl, (hsize_t)N_azimuth, (hsize_t)N_freq};
      hid_t   src_space   = H5Screate_simple(6, src_dims, NULL);

      hsize_t vds_start[6] = {0, 0, (hsize_t)slab_z_start[s], 0, 0, 0};
      hsize_t vds_count[6] = {(hsize_t)N_x,    (hsize_t)N_y,       (hsize_t)slab_nz[s],
                              (hsize_t)N_incl, (hsize_t)N_azimuth, (hsize_t)N_freq};
      H5Sselect_hyperslab(vspace, H5S_SELECT_SET, vds_start, NULL, vds_count, NULL);
      H5Pset_virtual(dcpl, vspace, src, names6[di], src_space);
      H5Sclose(src_space);
    }

    hid_t vspace_full = H5Screate_simple(6, vdims6, NULL);
    H5Dcreate2(grp, names6[di], f32, vspace_full, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Sclose(vspace_full);
    H5Pclose(dcpl);
    H5Sclose(vspace);
  }

  if (normalized_output) {
    for (int di = 0; di < 4; di++) {
      hid_t vspace = H5Screate_simple(5, vdims5, NULL);
      hid_t dcpl   = H5Pcreate(H5P_DATASET_CREATE);

      for (int s = 0; s < n_slabs; s++) {
        char src[PATH_MAX];
        build_zslab_filename(src, sizeof(src), dir, slab_z_start[s]);

        hsize_t src_dims[5] = {(hsize_t)N_x, (hsize_t)N_y, (hsize_t)slab_nz[s], (hsize_t)N_incl, (hsize_t)N_azimuth};
        hid_t   src_space   = H5Screate_simple(5, src_dims, NULL);

        hsize_t vds_start[5] = {0, 0, (hsize_t)slab_z_start[s], 0, 0};
        hsize_t vds_count[5] = {(hsize_t)N_x, (hsize_t)N_y, (hsize_t)slab_nz[s], (hsize_t)N_incl, (hsize_t)N_azimuth};
        H5Sselect_hyperslab(vspace, H5S_SELECT_SET, vds_start, NULL, vds_count, NULL);
        H5Pset_virtual(dcpl, vspace, src, names5[di], src_space);
        H5Sclose(src_space);
      }

      hid_t vspace_full = H5Screate_simple(5, vdims5, NULL);
      H5Dcreate2(grp, names5[di], f64, vspace_full, H5P_DEFAULT, dcpl, H5P_DEFAULT);
      H5Sclose(vspace_full);
      H5Pclose(dcpl);
      H5Sclose(vspace);
    }
  }

  H5Tclose(f32);
  H5Tclose(f64);
  H5Gclose(grp);
  H5Fclose(master);
}