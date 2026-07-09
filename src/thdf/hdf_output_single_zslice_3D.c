#include "hdf_output_single_zslice_3D.h"

#include <hdf5.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "hdf_fapl_from_env.h"
#include "output_utils.h"

int
THDF_create_zslab_file_single(MPI_Comm comm, const char *path, bool normalized_output, int N_x, int N_y, int N_z,
                              int N_incl, int N_azimuth, int N_frequencies) {
  hid_t file = -1, group = -1;
  hid_t fspace = -1, nspace = -1;
  hid_t f32 = -1, f64 = -1;
  hid_t dset = -1;
  int   rc   = 0;
  int   rank;
  MPI_Comm_rank(comm, &rank);

  hid_t fapl = build_fapl_from_env(comm);
  if (fapl < 0)
    return -1;

  file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  H5Pclose(fapl);
  if (file < 0) {
    if (rank == 0)
      fprintf(stderr, "[single] Cannot create: %s\n", path);
    return -1;
  }

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

  if (rank == 0)
    printf("[single] layout: contiguous (no chunking)\n");

  {
    hsize_t dims6[6] = {(hsize_t)N_x,    (hsize_t)N_y,       (hsize_t)N_z,
                        (hsize_t)N_incl, (hsize_t)N_azimuth, (hsize_t)N_frequencies};
    fspace           = H5Screate_simple(6, dims6, NULL);
    if (fspace < 0) {
      rc = -1;
      goto cleanup;
    }
    const char *names6[4] = {"I", "QI_pc", "UI_pc", "VI_pc"};
    for (int i = 0; i < 4; i++) {
      dset = H5Dcreate2(group, names6[i], f32, fspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (dset < 0) {
        if (rank == 0)
          fprintf(stderr, "[single] Error dataset %s\n", names6[i]);
        rc = -1;
        goto cleanup;
      }
      H5Dclose(dset);
      dset = -1;
    }
    H5Sclose(fspace);
    fspace = -1;
  }

  if (!normalized_output)
    goto sync;

  f64 = H5Tcopy(THDF_get_hdf_float_datatype());
  if (f64 < 0) {
    rc = -1;
    goto cleanup;
  }

  {
    hsize_t dims5[5] = {(hsize_t)N_x, (hsize_t)N_y, (hsize_t)N_z, (hsize_t)N_incl, (hsize_t)N_azimuth};
    nspace           = H5Screate_simple(5, dims5, NULL);
    if (nspace < 0) {
      rc = -1;
      goto cleanup;
    }
    const char *names5[4] = {"norm_multiplier_I", "norm_multiplier_QI_pc", "norm_multiplier_UI_pc",
                             "norm_multiplier_VI_pc"};
    for (int i = 0; i < 4; i++) {
      dset = H5Dcreate2(group, names5[i], f64, nspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (dset < 0) {
        if (rank == 0)
          fprintf(stderr, "[single] Error dataset %s\n", names5[i]);
        rc = -1;
        goto cleanup;
      }
      H5Dclose(dset);
      dset = -1;
    }
    H5Sclose(nspace);
    nspace = -1;
  }

sync:
  if (H5Fflush(file, H5F_SCOPE_GLOBAL) < 0) {
    if (rank == 0)
      fprintf(stderr, "[single] H5Fflush failed\n");
    rc = -1;
  }

cleanup:
  if (dset >= 0)
    H5Dclose(dset);
  if (fspace >= 0)
    H5Sclose(fspace);
  if (nspace >= 0)
    H5Sclose(nspace);
  if (f32 >= 0)
    H5Tclose(f32);
  if (f64 >= 0)
    H5Tclose(f64);
  if (group >= 0)
    H5Gclose(group);
  if (file >= 0)
    H5Fclose(file);
  return rc;
}

THDF_zslab_handler_t *
THDF_open_zslab_handler_single(hid_t file, bool normalized_output, int N_x, int N_y, int N_z_global, int N_incl,
                               int N_azimuth, int N_frequencies) {
  THDF_zslab_handler_t *h = calloc(1, sizeof(*h));
  if (!h)
    return NULL;

  h->file_id      = file;
  h->group_id     = -1;
  h->dataset_id_I = h->dataset_id_Q = h->dataset_id_U = h->dataset_id_V = -1;
  h->dataset_norm_multiplier_I = h->dataset_norm_multiplier_QI = -1;
  h->dataset_norm_multiplier_UI = h->dataset_norm_multiplier_VI = -1;
  h->datatype_f32_id = h->datatype_f64_id = -1;
  h->is_open                              = false;
  h->normalized_output                    = normalized_output;
  h->N_x                                  = N_x;
  h->N_y                                  = N_y;
  h->N_z                                  = N_z_global;
  h->N_incl                               = N_incl;
  h->N_azimuth                            = N_azimuth;
  h->N_frequencies                        = N_frequencies;

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
  if (h->dataset_id_I >= 0)
    H5Dclose(h->dataset_id_I);
  if (h->dataset_id_Q >= 0)
    H5Dclose(h->dataset_id_Q);
  if (h->dataset_id_U >= 0)
    H5Dclose(h->dataset_id_U);
  if (h->dataset_id_V >= 0)
    H5Dclose(h->dataset_id_V);
  if (h->dataset_norm_multiplier_I >= 0)
    H5Dclose(h->dataset_norm_multiplier_I);
  if (h->dataset_norm_multiplier_QI >= 0)
    H5Dclose(h->dataset_norm_multiplier_QI);
  if (h->dataset_norm_multiplier_UI >= 0)
    H5Dclose(h->dataset_norm_multiplier_UI);
  if (h->dataset_norm_multiplier_VI >= 0)
    H5Dclose(h->dataset_norm_multiplier_VI);
  if (h->datatype_f32_id >= 0)
    H5Tclose(h->datatype_f32_id);
  if (h->datatype_f64_id >= 0)
    H5Tclose(h->datatype_f64_id);
  if (h->group_id >= 0)
    H5Gclose(h->group_id);
  free(h);
  return NULL;
}

static int
write_stokes_hyperslab_single(hid_t dataset_id, hid_t datatype_id, hid_t memspace, hid_t dxpl, const hsize_t *start,
                              const hsize_t *count, int current_nz, const void *data, const char *name) {
  hid_t filespace = H5Dget_space(dataset_id);
  if (filespace < 0) {
    fprintf(stderr, "[single] H5Dget_space failed %s\n", name);
    return -1;
  }

  if (current_nz > 0) {
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);
  } else {
    H5Sselect_none(filespace);
    H5Sselect_none(memspace);
  }

  herr_t ret = H5Dwrite(dataset_id, datatype_id, memspace, filespace, dxpl, data);
  H5Sclose(filespace);
  if (ret < 0) {
    fprintf(stderr, "[single] H5Dwrite failed %s\n", name);
    return -1;
  }
  return 0;
}

int
THDF_write_zslab_block_single(THDF_zslab_handler_t *handler, const THDF_3D_field_t *field, int local_start_x,
                              int local_start_y, int global_start_z, int N_local_x, int N_local_y, int current_nz,
                              int N_incl, int N_azimuth, int N_frequencies) {
  if (!handler || !handler->is_open)
    return -1;

  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  if (dxpl < 0)
    return -1;
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

  hsize_t start[6] = {(hsize_t)local_start_x, (hsize_t)local_start_y, (hsize_t)global_start_z, 0, 0, 0};
  hsize_t count[6] = {(hsize_t)N_local_x, (hsize_t)N_local_y, (hsize_t)current_nz,
                      (hsize_t)N_incl,    (hsize_t)N_azimuth, (hsize_t)N_frequencies};

  hid_t memspace = H5Screate_simple(6, count, NULL);
  if (memspace < 0) {
    H5Pclose(dxpl);
    return -1;
  }

  const void *I_ptr  = (field && current_nz > 0) ? (const void *)field->stokes_I : NULL;
  const void *QI_ptr = (field && current_nz > 0) ? (const void *)field->stokes_QI : NULL;
  const void *UI_ptr = (field && current_nz > 0) ? (const void *)field->stokes_UI : NULL;
  const void *VI_ptr = (field && current_nz > 0) ? (const void *)field->stokes_VI : NULL;

  int rc = 0;
  rc |= write_stokes_hyperslab_single(handler->dataset_id_I, handler->datatype_f32_id, memspace, dxpl, start, count,
                                      current_nz, I_ptr, "I");
  rc |= write_stokes_hyperslab_single(handler->dataset_id_Q, handler->datatype_f32_id, memspace, dxpl, start, count,
                                      current_nz, QI_ptr, "QI");
  rc |= write_stokes_hyperslab_single(handler->dataset_id_U, handler->datatype_f32_id, memspace, dxpl, start, count,
                                      current_nz, UI_ptr, "UI");
  rc |= write_stokes_hyperslab_single(handler->dataset_id_V, handler->datatype_f32_id, memspace, dxpl, start, count,
                                      current_nz, VI_ptr, "VI");

#ifndef NDEBUG
  {
    uint32_t lc, gc;
    H5Pget_mpio_no_collective_cause(dxpl, &lc, &gc);
    if (gc) {
      int r;
      MPI_Comm_rank(MPI_COMM_WORLD, &r);
      if (!r)
        fprintf(stderr, "[single] collective fallback 0x%x\n", gc);
    }
  }
#endif

  if (!handler->normalized_output) {
    H5Sclose(memspace);
    H5Pclose(dxpl);
    return rc;
  }

  hsize_t nstart[5] = {(hsize_t)local_start_x, (hsize_t)local_start_y, (hsize_t)global_start_z, 0, 0};
  hsize_t ncount[5] = {(hsize_t)N_local_x, (hsize_t)N_local_y, (hsize_t)current_nz, (hsize_t)N_incl,
                       (hsize_t)N_azimuth};
  hid_t   nmemspace = H5Screate_simple(5, ncount, NULL);
  if (nmemspace < 0) {
    H5Sclose(memspace);
    H5Pclose(dxpl);
    return -1;
  }

  const void *nI  = (field && current_nz > 0) ? (const void *)field->norm_multiplier_I : NULL;
  const void *nQI = (field && current_nz > 0) ? (const void *)field->norm_multiplier_QI : NULL;
  const void *nUI = (field && current_nz > 0) ? (const void *)field->norm_multiplier_UI : NULL;
  const void *nVI = (field && current_nz > 0) ? (const void *)field->norm_multiplier_VI : NULL;

  rc |= write_stokes_hyperslab_single(handler->dataset_norm_multiplier_I, handler->datatype_f64_id, nmemspace, dxpl,
                                      nstart, ncount, current_nz, nI, "norm_I");
  rc |= write_stokes_hyperslab_single(handler->dataset_norm_multiplier_QI, handler->datatype_f64_id, nmemspace, dxpl,
                                      nstart, ncount, current_nz, nQI, "norm_QI");
  rc |= write_stokes_hyperslab_single(handler->dataset_norm_multiplier_UI, handler->datatype_f64_id, nmemspace, dxpl,
                                      nstart, ncount, current_nz, nUI, "norm_UI");
  rc |= write_stokes_hyperslab_single(handler->dataset_norm_multiplier_VI, handler->datatype_f64_id, nmemspace, dxpl,
                                      nstart, ncount, current_nz, nVI, "norm_VI");

  H5Sclose(nmemspace);
  H5Sclose(memspace);
  H5Pclose(dxpl);
  return rc;
}

void
THDF_close_zslab_handler_single(THDF_zslab_handler_t *h) {
  if (!h || !h->is_open)
    return;
  if (h->dataset_id_I >= 0)
    H5Dclose(h->dataset_id_I);
  if (h->dataset_id_Q >= 0)
    H5Dclose(h->dataset_id_Q);
  if (h->dataset_id_U >= 0)
    H5Dclose(h->dataset_id_U);
  if (h->dataset_id_V >= 0)
    H5Dclose(h->dataset_id_V);
  if (h->dataset_norm_multiplier_I >= 0)
    H5Dclose(h->dataset_norm_multiplier_I);
  if (h->dataset_norm_multiplier_QI >= 0)
    H5Dclose(h->dataset_norm_multiplier_QI);
  if (h->dataset_norm_multiplier_UI >= 0)
    H5Dclose(h->dataset_norm_multiplier_UI);
  if (h->dataset_norm_multiplier_VI >= 0)
    H5Dclose(h->dataset_norm_multiplier_VI);
  if (h->datatype_f32_id >= 0)
    H5Tclose(h->datatype_f32_id);
  if (h->datatype_f64_id >= 0)
    H5Tclose(h->datatype_f64_id);
  if (h->group_id >= 0)
    H5Gclose(h->group_id);
  h->is_open = false;
  free(h);
}