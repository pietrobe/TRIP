#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf_atmos_continuum.h"

////////////////////////////////////////////////////////////////////////
// Continuum dataset functions
////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
// create_hdf_continuum_compound_type
/////////////////////////////////////////////////////////////////////////
static int
write_or_update_int_attribute(hid_t object_id, const char *name, int value) {
  hid_t  attr_id  = H5I_INVALID_HID;
  hid_t  space_id = H5I_INVALID_HID;
  herr_t status;

  if (H5Aexists(object_id, name) > 0) {
    attr_id = H5Aopen(object_id, name, H5P_DEFAULT);
  } else {
    hsize_t dims[1] = {1};
    space_id        = H5Screate_simple(1, dims, NULL);
    if (space_id < 0)
      return -1;

    attr_id = H5Acreate2(object_id, name, H5T_NATIVE_INT, space_id, H5P_DEFAULT, H5P_DEFAULT);
  }

  if (attr_id < 0) {
    if (space_id >= 0)
      H5Sclose(space_id);
    return -1;
  }

  status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);

  H5Aclose(attr_id);
  if (space_id >= 0)
    H5Sclose(space_id);

  return (status < 0) ? -1 : 0;
}

static int
create_or_open_u16_4d_dataset(hid_t group_id, const char *dataset_name, const hsize_t dims[4]) {
  hid_t dspace     = H5I_INVALID_HID;
  hid_t dset       = H5I_INVALID_HID;
  hid_t dset_space = H5I_INVALID_HID;
  hid_t dset_type  = H5I_INVALID_HID;

  dspace = H5Screate_simple(4, dims, NULL);
  if (dspace < 0) {
    fprintf(stderr, "Error creating dataspace for %s\n", dataset_name);
    return -1;
  }

  H5E_BEGIN_TRY {
    dset = H5Dcreate2(group_id, dataset_name, H5T_NATIVE_UINT16, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  H5E_END_TRY;

  if (dset < 0) {
    H5E_BEGIN_TRY { dset = H5Dopen2(group_id, dataset_name, H5P_DEFAULT); }
    H5E_END_TRY;
    if (dset < 0) {
      fprintf(stderr, "Error creating/opening dataset %s\n", dataset_name);
      H5Sclose(dspace);
      return -1;
    }

    dset_space = H5Dget_space(dset);
    if (dset_space < 0) {
      fprintf(stderr, "Error reading dataspace for %s\n", dataset_name);
      H5Dclose(dset);
      H5Sclose(dspace);
      return -1;
    }

    int existing_ndims = H5Sget_simple_extent_ndims(dset_space);
    if (existing_ndims != 4) {
      fprintf(stderr, "Dataset %s is not 4D\n", dataset_name);
      H5Sclose(dset_space);
      H5Dclose(dset);
      H5Sclose(dspace);
      return -1;
    }

    hsize_t existing_dims[4] = {0, 0, 0, 0};
    H5Sget_simple_extent_dims(dset_space, existing_dims, NULL);
    for (int axis = 0; axis < 4; axis++) {
      if (existing_dims[axis] != dims[axis]) {
        fprintf(stderr, "Dataset %s has incompatible dimensions\n", dataset_name);
        H5Sclose(dset_space);
        H5Dclose(dset);
        H5Sclose(dspace);
        return -1;
      }
    }

    dset_type = H5Dget_type(dset);
    if (dset_type < 0 || H5Tget_class(dset_type) != H5T_INTEGER || H5Tget_sign(dset_type) != H5T_SGN_NONE ||
        H5Tget_size(dset_type) != sizeof(uint16_t)) {
      fprintf(stderr, "Dataset %s is not uint16\n", dataset_name);
      if (dset_type >= 0)
        H5Tclose(dset_type);
      H5Sclose(dset_space);
      H5Dclose(dset);
      H5Sclose(dspace);
      return -1;
    }

    H5Tclose(dset_type);
    H5Sclose(dset_space);
  }

  H5Dclose(dset);
  H5Sclose(dspace);
  return 0;
}

static int
read_u16_4d_dataset_slice(hid_t group_id, const char *dataset_name, const hsize_t start[4], const hsize_t count[4],
                          const hsize_t mem_dims[4], uint16_t *buffer) {
  hid_t  dset      = H5I_INVALID_HID;
  hid_t  dspace    = H5I_INVALID_HID;
  hid_t  mspace    = H5I_INVALID_HID;
  herr_t read_stat = -1;

  H5E_BEGIN_TRY { dset = H5Dopen2(group_id, dataset_name, H5P_DEFAULT); }
  H5E_END_TRY;
  if (dset < 0) {
    fprintf(stderr, "Error opening dataset %s\n", dataset_name);
    return -1;
  }

  dspace = H5Dget_space(dset);
  if (dspace < 0) {
    fprintf(stderr, "Error getting dataspace for %s\n", dataset_name);
    H5Dclose(dset);
    return -1;
  }

  if (H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
    fprintf(stderr, "Error selecting hyperslab for %s\n", dataset_name);
    H5Sclose(dspace);
    H5Dclose(dset);
    return -1;
  }

  mspace = H5Screate_simple(4, mem_dims, NULL);
  if (mspace < 0) {
    fprintf(stderr, "Error creating memory dataspace for %s\n", dataset_name);
    H5Sclose(dspace);
    H5Dclose(dset);
    return -1;
  }

  read_stat = H5Dread(dset, H5T_NATIVE_UINT16, mspace, dspace, H5P_DEFAULT, buffer);

  H5Sclose(mspace);
  H5Sclose(dspace);
  H5Dclose(dset);

  if (read_stat < 0) {
    fprintf(stderr, "Error reading dataset %s\n", dataset_name);
    return -1;
  }

  return 0;
}

/////////////////////////////////////////////////////////////////////////
// THDF_create_continuum_4D_matrices
/////////////////////////////////////////////////////////////////////////
int
THDF_create_continuum_4D_matrices(hid_t file, int N_x, int N_y, int N_z, int N_frequencies) {
  hid_t cont_group = H5I_INVALID_HID;

  if (N_x <= 0 || N_y <= 0 || N_z <= 0 || N_frequencies <= 0) {
    fprintf(stderr, "Invalid continuum 4D dimensions\n");
    return -1;
  }

  H5E_BEGIN_TRY { cont_group = H5Gopen2(file, "/continuum", H5P_DEFAULT); }
  H5E_END_TRY;

  if (cont_group < 0) {
    cont_group = H5Gcreate2(file, "/continuum", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (cont_group < 0) {
      fprintf(stderr, "Error creating continuum group\n");
      return -1;
    }
  }

  if (write_or_update_int_attribute(cont_group, "N_x", N_x) < 0 ||
      write_or_update_int_attribute(cont_group, "N_y", N_y) < 0 ||
      write_or_update_int_attribute(cont_group, "N_z", N_z) < 0 ||
      write_or_update_int_attribute(cont_group, "N_frequencies", N_frequencies) < 0) {
    fprintf(stderr, "Error writing continuum metadata\n");
    H5Gclose(cont_group);
    return -1;
  }

  const hsize_t dims[4] = {(hsize_t)N_x, (hsize_t)N_y, (hsize_t)N_z, (hsize_t)N_frequencies};

  if (create_or_open_u16_4d_dataset(cont_group, "c_scat_opacity_sigma_c", dims) < 0 ||
      create_or_open_u16_4d_dataset(cont_group, "c_therm_opacity_k_c", dims) < 0 ||
      create_or_open_u16_4d_dataset(cont_group, "c_therm_emissivity_epsilon_c", dims) < 0) {
    H5Gclose(cont_group);
    return -1;
  }

  H5Gclose(cont_group);
  return 0;
}

/////////////////////////////////////////////////////////////////////////
// create_hdf_continuum_compound_type
/////////////////////////////////////////////////////////////////////////
hid_t
create_hdf_continuum_compound_type(const int N_frequencies) {
  hid_t continuum_type = H5Tcreate(H5T_COMPOUND, sizeof(THDF_continuum_t));
  if (continuum_type < 0)
    return H5I_INVALID_HID;

  /* Scalar fields -------------------------------------------------------- */
  H5Tinsert(continuum_type, "index_i", HOFFSET(THDF_continuum_t, index_i), H5T_NATIVE_INT);
  H5Tinsert(continuum_type, "index_j", HOFFSET(THDF_continuum_t, index_j), H5T_NATIVE_INT);
  H5Tinsert(continuum_type, "index_k", HOFFSET(THDF_continuum_t, index_k), H5T_NATIVE_INT);

  H5Tinsert(continuum_type, "min_c_scat_opacity_sigma_c", HOFFSET(THDF_continuum_t, min_c_scat_opacity_sigma_c),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(continuum_type, "max_Vm_c_scat_opacity_sigma_c", HOFFSET(THDF_continuum_t, max_Vm_c_scat_opacity_sigma_c),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(continuum_type, "min_c_therm_opacity_k_c", HOFFSET(THDF_continuum_t, min_c_therm_opacity_k_c),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(continuum_type, "max_Vm_c_therm_opacity_k_c", HOFFSET(THDF_continuum_t, max_Vm_c_therm_opacity_k_c),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(continuum_type, "min_c_therm_emissivity_epsilon_c",
            HOFFSET(THDF_continuum_t, min_c_therm_emissivity_epsilon_c), H5T_NATIVE_DOUBLE);
  H5Tinsert(continuum_type, "max_Vm_c_therm_emissivity_epsilon_c",
            HOFFSET(THDF_continuum_t, max_Vm_c_therm_emissivity_epsilon_c), H5T_NATIVE_DOUBLE);

  (void)N_frequencies;
  /*
   * Pointer fields cannot be represented as inline fixed-size arrays in an
   * HDF5 compound type. Keep normalization scalars in this committed type and
   * store spectral arrays in dedicated datasets.
   */
  return continuum_type;
}

/////////////////////////////////////////////////////////////////////////
// make_continuum_hdf5_mpi
/////////////////////////////////////////////////////////////////////////
int
make_continuum_hdf5_mpi(hid_t file,                                      //
                        int N_x, int N_y, int N_z, int N_frequencies) {  //

  hid_t cont_group     = H5I_INVALID_HID;
  hid_t continuum_type = H5I_INVALID_HID;
  hid_t committed_type = H5I_INVALID_HID;
  hid_t dspace         = H5I_INVALID_HID;
  hid_t dset           = H5I_INVALID_HID;
  hid_t dset_type      = H5I_INVALID_HID;
  // herr_t commit_status;

  if (N_x <= 0 || N_y <= 0 || N_z <= 0 || N_frequencies <= 0) {
    fprintf(stderr, "Invalid continuum dimensions\n");
    return -1;
  }

  H5E_BEGIN_TRY { cont_group = H5Gopen2(file, "/continuum", H5P_DEFAULT); }
  H5E_END_TRY;

  if (cont_group < 0) {
    cont_group = H5Gcreate2(file, "/continuum", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (cont_group < 0) {
      fprintf(stderr, "Error creating continuum group\n");
      return -1;
    }
  }

  continuum_type = create_hdf_continuum_compound_type(N_frequencies);
  if (continuum_type < 0) {
    fprintf(stderr, "Error creating continuum datatype\n");
    H5Gclose(cont_group);
    return -1;
  }

  if (write_or_update_int_attribute(cont_group, "N_x", N_x) < 0 ||
      write_or_update_int_attribute(cont_group, "N_y", N_y) < 0 ||
      write_or_update_int_attribute(cont_group, "N_z", N_z) < 0 ||
      write_or_update_int_attribute(cont_group, "N_frequencies", N_frequencies) < 0) {
    fprintf(stderr, "Error writing continuum metadata\n");
    if (committed_type >= 0)
      H5Tclose(committed_type);
    H5Tclose(continuum_type);
    H5Gclose(cont_group);
    return -1;
  }

  if (THDF_create_continuum_4D_matrices(file, N_x, N_y, N_z, N_frequencies) < 0) {
    fprintf(stderr, "Error creating 4D continuum matrices\n");
    if (committed_type >= 0)
      H5Tclose(committed_type);
    H5Tclose(continuum_type);
    H5Gclose(cont_group);
    return -1;
  }

  {
    hsize_t dims[3] = {(hsize_t)N_x, (hsize_t)N_y, (hsize_t)N_z};
    dspace          = H5Screate_simple(3, dims, NULL);
    if (dspace < 0) {
      fprintf(stderr, "Error creating continuum dataspace\n");
      if (committed_type >= 0)
        H5Tclose(committed_type);
      H5Tclose(continuum_type);
      H5Gclose(cont_group);
      return -1;
    }

    H5E_BEGIN_TRY {
      dset = H5Dcreate2(cont_group, "continuum_norm", continuum_type, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    H5E_END_TRY;

    if (dset < 0) {
      H5E_BEGIN_TRY { dset = H5Dopen2(cont_group, "continuum_norm", H5P_DEFAULT); }
      H5E_END_TRY;
      if (dset < 0) {
        fprintf(stderr, "Error creating/opening continuum dataset\n");
        H5Sclose(dspace);
        if (committed_type >= 0)
          H5Tclose(committed_type);
        H5Tclose(continuum_type);
        H5Gclose(cont_group);
        return -1;
      }

      dset_type = H5Dget_type(dset);
      if (dset_type < 0 || H5Tequal(dset_type, continuum_type) <= 0) {
        fprintf(stderr, "Existing continuum dataset has incompatible datatype\n");
        if (dset_type >= 0)
          H5Tclose(dset_type);
        H5Dclose(dset);
        H5Sclose(dspace);
        if (committed_type >= 0)
          H5Tclose(committed_type);
        H5Tclose(continuum_type);
        H5Gclose(cont_group);
        return -1;
      }
      H5Tclose(dset_type);
    }

    H5Dclose(dset);
    H5Sclose(dspace);
  }

  if (committed_type >= 0)
    H5Tclose(committed_type);
  H5Tclose(continuum_type);
  H5Gclose(cont_group);

  return 0;
}

////////////////////////////////////////////////////////////////////////////
// normalize_continuum_data
////////////////////////////////////////////////////////////////////////////
THDF_continuum_t
normalize_continuum_data(const double *sigma_c, const double *epsilon_c, const double *c_therm_emissivity_epsilon_c,
                         int N_frequencies, uint16_t *norm_sigma_c, uint16_t *norm_epsilon_c,
                         uint16_t *norm_c_therm_emissivity_epsilon_c) {

  double min_sigma_c                        = sigma_c[0];
  double max_sigma_c                        = sigma_c[0];
  double min_epsilon_c                      = epsilon_c[0];
  double max_epsilon_c                      = epsilon_c[0];
  double min_c_therm_emissivity_epsilon_c = c_therm_emissivity_epsilon_c[0];
  double max_c_therm_emissivity_epsilon_c = c_therm_emissivity_epsilon_c[0];

  for (int i = 0; i < N_frequencies; i++) {
    min_sigma_c                        = sigma_c[i] < min_sigma_c ? sigma_c[i] : min_sigma_c;
    max_sigma_c                        = sigma_c[i] > max_sigma_c ? sigma_c[i] : max_sigma_c;
    min_epsilon_c                      = epsilon_c[i] < min_epsilon_c ? epsilon_c[i] : min_epsilon_c;
    max_epsilon_c                      = epsilon_c[i] > max_epsilon_c ? epsilon_c[i] : max_epsilon_c;
    min_c_therm_emissivity_epsilon_c = c_therm_emissivity_epsilon_c[i] < min_c_therm_emissivity_epsilon_c
                                             ? c_therm_emissivity_epsilon_c[i]
                                             : min_c_therm_emissivity_epsilon_c;
    max_c_therm_emissivity_epsilon_c = c_therm_emissivity_epsilon_c[i] > max_c_therm_emissivity_epsilon_c
                                             ? c_therm_emissivity_epsilon_c[i]
                                             : max_c_therm_emissivity_epsilon_c;
  }

  THDF_continuum_t cont;
  cont.min_c_scat_opacity_sigma_c            = min_sigma_c;
  cont.max_Vm_c_scat_opacity_sigma_c         = max_sigma_c - min_sigma_c;
  cont.min_c_therm_opacity_k_c               = min_epsilon_c;
  cont.max_Vm_c_therm_opacity_k_c            = max_epsilon_c - min_epsilon_c;
  cont.min_c_therm_emissivity_epsilon_c    = min_c_therm_emissivity_epsilon_c;
  cont.max_Vm_c_therm_emissivity_epsilon_c = max_c_therm_emissivity_epsilon_c - min_c_therm_emissivity_epsilon_c;

  for (int i = 0; i < N_frequencies; i++) {

    const double normed_sigma_c   = (sigma_c[i] - min_sigma_c) / cont.max_Vm_c_scat_opacity_sigma_c;
    const double normed_epsilon_c = (epsilon_c[i] - min_epsilon_c) / cont.max_Vm_c_therm_opacity_k_c;
    const double normed_c_therm_emissivity_epsilon_c =
        (c_therm_emissivity_epsilon_c[i] - min_c_therm_emissivity_epsilon_c) /
        cont.max_Vm_c_therm_emissivity_epsilon_c;

    norm_sigma_c[i]   = (uint16_t)((normed_sigma_c < 0.0 ? 0.0 : normed_sigma_c) * UINT16_MAX);
    norm_epsilon_c[i] = (uint16_t)((normed_epsilon_c < 0.0 ? 0.0 : normed_epsilon_c) * UINT16_MAX);
    norm_c_therm_emissivity_epsilon_c[i] =
        (uint16_t)((normed_c_therm_emissivity_epsilon_c < 0.0 ? 0.0 : normed_c_therm_emissivity_epsilon_c) *
                   UINT16_MAX);
  }

  return cont;
}

/////////////////////////////////////////////////////////////////////////
// denormalize_continuum_data
/////////////////////////////////////////////////////////////////////////
int
denormalize_continuum_data(const THDF_continuum_t cont, const uint16_t *norm_sigma_c, const uint16_t *norm_epsilon_c,
                           const uint16_t *norm_c_therm_emissivity_epsilon_c, int N_frequencies, double *sigma_c,
                           double *epsilon_c, double *c_therm_emissivity_epsilon_c) {

  double inv_UINT16_MAX = 1.0 / (double)UINT16_MAX;

  for (int i = 0; i < N_frequencies; i++) {
    sigma_c[i] = cont.min_c_scat_opacity_sigma_c +             //
                 ((double)norm_sigma_c[i] * inv_UINT16_MAX) *  //
                     cont.max_Vm_c_scat_opacity_sigma_c;       //

    epsilon_c[i] = cont.min_c_therm_opacity_k_c +                  //
                   ((double)norm_epsilon_c[i] * inv_UINT16_MAX) *  //
                       cont.max_Vm_c_therm_opacity_k_c;            //

    c_therm_emissivity_epsilon_c[i] = cont.min_c_therm_emissivity_epsilon_c +                            //
                                        ((double)norm_c_therm_emissivity_epsilon_c[i] * inv_UINT16_MAX) *  //
                                            cont.max_Vm_c_therm_emissivity_epsilon_c;
  }

  return 0;
}

/////////////////////////////////////////////////////////////////////////
// write_normalized_continuum_to_hdf5
/////////////////////////////////////////////////////////////////////////
int
write_normalized_continuum_to_hdf5(hid_t     file,                                           //
                                   const int N_frequencies,                                  //
                                   const int ix, const int iy, const int iz,                 //
                                   const int count_x, const int count_y, const int count_z,  //
                                   const THDF_continuum_t *normalized_cont,                  //
                                   uint16_t *norm_sigma_c, uint16_t *norm_epsilon_c,
                                   uint16_t *norm_c_therm_emissivity_epsilon_c) {  //

  /* Write norm_sigma_c, norm_epsilon_c, norm_c_therm_emissivity_epsilon_c to the appropriate slice of the 4D datasets
   */

  /* Open continuum group */

  int exit_code = EXIT_SUCCESS;

  hid_t cont_group = H5I_INVALID_HID;
  H5E_BEGIN_TRY { cont_group = H5Gopen2(file, "/continuum", H5P_DEFAULT); }
  H5E_END_TRY;
  if (cont_group < 0) {
    fprintf(stderr, "Error opening continuum group\n");
    return -1;
  }

  /* Setup hyperslab selection for 4D datasets: full local block in x,y,z,f */
  hsize_t start[4] = {(hsize_t)ix, (hsize_t)iy, (hsize_t)iz, 0};
  hsize_t count[4] = {(hsize_t)count_x, (hsize_t)count_y, (hsize_t)count_z,
                      (hsize_t)N_frequencies}; /* Placeholder; will use actual dimension */

  hid_t             dspace       = H5I_INVALID_HID;
  hid_t             mspace       = H5I_INVALID_HID;
  hid_t             dset         = H5I_INVALID_HID;
  hid_t             mem_type     = H5I_INVALID_HID;
  THDF_continuum_t *cont_block   = NULL;
  herr_t            write_status = 0;

  /* Get actual N_frequencies from file attributes */
  int   N_frequencies_actual = 0;
  hid_t attr_id              = H5I_INVALID_HID;
  H5E_BEGIN_TRY {
    attr_id = H5Aopen(cont_group, "N_frequencies", H5P_DEFAULT);
    if (attr_id >= 0) {
      H5Aread(attr_id, H5T_NATIVE_INT, &N_frequencies_actual);
      H5Aclose(attr_id);
    }
  }
  H5E_END_TRY;

  if (N_frequencies_actual <= 0) {
    fprintf(stderr, "Error reading N_frequencies attribute\n");
    H5Gclose(cont_group);
    return -1;
  }

  count[3] = (hsize_t)N_frequencies_actual;

  hsize_t mem_dims_4d[4] = {(hsize_t)count_x, (hsize_t)count_y, (hsize_t)count_z, (hsize_t)N_frequencies_actual};

  /* Write c_scat_opacity_sigma_c dataset */
  H5E_BEGIN_TRY { dset = H5Dopen2(cont_group, "c_scat_opacity_sigma_c", H5P_DEFAULT); }
  H5E_END_TRY;
  if (dset >= 0) {
    dspace = H5Dget_space(dset);
    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, NULL, count, NULL);
    mspace       = H5Screate_simple(4, mem_dims_4d, NULL);
    write_status = H5Dwrite(dset, H5T_NATIVE_UINT16, mspace, dspace, H5P_DEFAULT, norm_sigma_c);
    H5Sclose(mspace);
    H5Sclose(dspace);
    H5Dclose(dset);
    if (write_status < 0) {
      fprintf(stderr, "Error writing c_scat_opacity_sigma_c\n");
      H5Gclose(cont_group);
      return -1;
    }
  }

  /* Write c_therm_opacity_k_c dataset */
  H5E_BEGIN_TRY { dset = H5Dopen2(cont_group, "c_therm_opacity_k_c", H5P_DEFAULT); }
  H5E_END_TRY;
  if (dset >= 0) {
    dspace = H5Dget_space(dset);
    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, NULL, count, NULL);
    mspace       = H5Screate_simple(4, mem_dims_4d, NULL);
    write_status = H5Dwrite(dset, H5T_NATIVE_UINT16, mspace, dspace, H5P_DEFAULT, norm_epsilon_c);
    H5Sclose(mspace);
    H5Sclose(dspace);
    H5Dclose(dset);
    if (write_status < 0) {
      fprintf(stderr, "Error writing c_therm_opacity_k_c\n");
      H5Gclose(cont_group);
      return -1;
    }
  }

  /* Write c_therm_emissivity_epsilon_c dataset */
  H5E_BEGIN_TRY { dset = H5Dopen2(cont_group, "c_therm_emissivity_epsilon_c", H5P_DEFAULT); }
  H5E_END_TRY;
  if (dset >= 0) {
    dspace = H5Dget_space(dset);
    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, NULL, count, NULL);
    mspace       = H5Screate_simple(4, mem_dims_4d, NULL);
    write_status = H5Dwrite(dset, H5T_NATIVE_UINT16, mspace, dspace, H5P_DEFAULT, norm_c_therm_emissivity_epsilon_c);
    H5Sclose(mspace);
    H5Sclose(dspace);
    H5Dclose(dset);
    if (write_status < 0) {
      fprintf(stderr, "Error writing c_therm_emissivity_epsilon_c\n");
      H5Gclose(cont_group);
      return -1;
    }
  }

  /* Write cont_to_write to continuum_norm dataset at (ix, iy, iz) */
  H5E_BEGIN_TRY { dset = H5Dopen2(cont_group, "continuum_norm", H5P_DEFAULT); }
  H5E_END_TRY;
  if (dset < 0) {
    fprintf(stderr, "Error opening continuum_norm dataset\n");
    H5Gclose(cont_group);
    return -1;
  }

  {
    hsize_t      start_3d[3] = {(hsize_t)ix, (hsize_t)iy, (hsize_t)iz};
    hsize_t      count_3d[3] = {(hsize_t)count_x, (hsize_t)count_y, (hsize_t)count_z};
    hsize_t      mem_dims[3] = {(hsize_t)count_x, (hsize_t)count_y, (hsize_t)count_z};
    const size_t block_size  = (size_t)count_x * (size_t)count_y * (size_t)count_z;

    cont_block = (THDF_continuum_t *)malloc(block_size * sizeof(*cont_block));
    if (cont_block == NULL) {
      fprintf(stderr, "Error allocating continuum block buffer\n");
      H5Dclose(dset);
      H5Gclose(cont_group);
      return -1;
    }

    for (int dx = 0; dx < count_x; dx++) {
      for (int dy = 0; dy < count_y; dy++) {
        for (int dz = 0; dz < count_z; dz++) {
          const size_t local_idx =
              (size_t)dx * (size_t)count_y * (size_t)count_z + (size_t)dy * (size_t)count_z + (size_t)dz;
          cont_block[local_idx] = normalized_cont[local_idx];
          // cont_block[local_idx].index_i = ix + dx;
          // cont_block[local_idx].index_j = iy + dy;
          // cont_block[local_idx].index_k = iz + dz;
        }
      }
    }

    dspace = H5Dget_space(dset);
    if (dspace < 0) {
      fprintf(stderr, "Error getting dataspace for continuum_norm\n");
      free(cont_block);
      H5Dclose(dset);
      H5Gclose(cont_group);
      return -1;
    }

    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start_3d, NULL, count_3d, NULL);
    mspace = H5Screate_simple(3, mem_dims, NULL);
    if (mspace < 0) {
      fprintf(stderr, "Error creating memory dataspace for continuum_norm\n");
      free(cont_block);
      H5Sclose(dspace);
      H5Dclose(dset);
      H5Gclose(cont_group);
      return -1;
    }

    mem_type = create_hdf_continuum_compound_type(N_frequencies_actual);
    if (mem_type < 0) {
      fprintf(stderr, "Error creating memory datatype for continuum_norm\n");
      free(cont_block);
      H5Sclose(mspace);
      H5Sclose(dspace);
      H5Dclose(dset);
      H5Gclose(cont_group);
      return -1;
    }

    write_status = H5Dwrite(dset, mem_type, mspace, dspace, H5P_DEFAULT, cont_block);

    free(cont_block);
    H5Tclose(mem_type);
    H5Sclose(mspace);
    H5Sclose(dspace);
    H5Dclose(dset);

    if (write_status < 0) {
      fprintf(stderr, "Error writing continuum_norm\n");
      H5Gclose(cont_group);
      return -1;
    }
  }

  H5Gclose(cont_group);
  return exit_code;
}

int
read_continuum_from_hdf5(hid_t     file,                                                                //
                         const int N_frequencies,                                                       //
                         const int ix, const int iy, const int iz,                                      //
                         const int count_x, const int count_y, const int count_z,                       //
                         THDF_continuum_t *cont_to_fill,                                                //
                         double *sigma_c, double *epsilon_c, double *c_therm_emissivity_epsilon_c) {  //
                                                                                                        //
  if (N_frequencies <= 0 || count_x <= 0 || count_y <= 0 || count_z <= 0 || cont_to_fill == NULL || sigma_c == NULL ||
      epsilon_c == NULL || c_therm_emissivity_epsilon_c == NULL) {
    fprintf(stderr, "Invalid arguments for reading continuum data\n");
    return -1;
  }

  hid_t cont_group = H5I_INVALID_HID;
  hid_t dset       = H5I_INVALID_HID;

  H5E_BEGIN_TRY { cont_group = H5Gopen2(file, "/continuum", H5P_DEFAULT); }
  H5E_END_TRY;
  if (cont_group < 0) {
    fprintf(stderr, "Error opening continuum group\n");
    return -1;
  }

  H5E_BEGIN_TRY { dset = H5Dopen2(cont_group, "continuum_norm", H5P_DEFAULT); }
  H5E_END_TRY;
  if (dset < 0) {
    fprintf(stderr, "Error opening continuum_norm dataset\n");
    H5Gclose(cont_group);
    return -1;
  }

  {
    hsize_t      start_3d[3]                         = {(hsize_t)ix, (hsize_t)iy, (hsize_t)iz};
    hsize_t      count_3d[3]                         = {(hsize_t)count_x, (hsize_t)count_y, (hsize_t)count_z};
    hsize_t      mem_dims_3d[3]                      = {(hsize_t)count_x, (hsize_t)count_y, (hsize_t)count_z};
    hid_t        dspace                              = H5I_INVALID_HID;
    hid_t        mspace                              = H5I_INVALID_HID;
    hid_t        mem_type                            = H5I_INVALID_HID;
    const size_t block_size                          = (size_t)count_x * (size_t)count_y * (size_t)count_z;
    const size_t vec_size                            = block_size * (size_t)N_frequencies;
    uint16_t    *norm_sigma_c                        = NULL;
    uint16_t    *norm_epsilon_c                      = NULL;
    uint16_t    *norm_c_therm_emissivity_epsilon_c = NULL;

    dspace = H5Dget_space(dset);
    if (dspace < 0) {
      fprintf(stderr, "Error getting dataspace for continuum_norm\n");
      H5Dclose(dset);
      H5Gclose(cont_group);
      return -1;
    }

    if (H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start_3d, NULL, count_3d, NULL) < 0) {
      fprintf(stderr, "Error selecting hyperslab for continuum_norm\n");
      H5Sclose(dspace);
      H5Dclose(dset);
      H5Gclose(cont_group);
      return -1;
    }

    mspace = H5Screate_simple(3, mem_dims_3d, NULL);
    if (mspace < 0) {
      fprintf(stderr, "Error creating memory dataspace for continuum_norm\n");
      H5Sclose(dspace);
      H5Dclose(dset);
      H5Gclose(cont_group);
      return -1;
    }

    mem_type = create_hdf_continuum_compound_type(N_frequencies);
    if (mem_type < 0) {
      fprintf(stderr, "Error creating memory datatype for continuum_norm\n");
      H5Sclose(mspace);
      H5Sclose(dspace);
      H5Dclose(dset);
      H5Gclose(cont_group);
      return -1;
    }

    if (H5Dread(dset, mem_type, mspace, dspace, H5P_DEFAULT, cont_to_fill) < 0) {
      fprintf(stderr, "Error reading continuum_norm\n");
      H5Tclose(mem_type);
      H5Sclose(mspace);
      H5Sclose(dspace);
      H5Dclose(dset);
      H5Gclose(cont_group);
      return -1;
    }

    H5Tclose(mem_type);
    H5Sclose(mspace);
    H5Sclose(dspace);
    H5Dclose(dset);

    norm_sigma_c                        = (uint16_t *)malloc(vec_size * sizeof(uint16_t));
    norm_epsilon_c                      = (uint16_t *)malloc(vec_size * sizeof(uint16_t));
    norm_c_therm_emissivity_epsilon_c = (uint16_t *)malloc(vec_size * sizeof(uint16_t));

    if (norm_sigma_c == NULL || norm_epsilon_c == NULL || norm_c_therm_emissivity_epsilon_c == NULL) {
      fprintf(stderr, "Error allocating normalization buffers\n");
      free(norm_sigma_c);
      free(norm_epsilon_c);
      free(norm_c_therm_emissivity_epsilon_c);
      H5Gclose(cont_group);
      return -1;
    }

    {
      hsize_t start_4d[4]    = {(hsize_t)ix, (hsize_t)iy, (hsize_t)iz, 0};
      hsize_t count_4d[4]    = {(hsize_t)count_x, (hsize_t)count_y, (hsize_t)count_z, (hsize_t)N_frequencies};
      hsize_t mem_dims_4d[4] = {(hsize_t)count_x, (hsize_t)count_y, (hsize_t)count_z, (hsize_t)N_frequencies};

      const int f1 = read_u16_4d_dataset_slice(cont_group,  //
                                               "c_scat_opacity_sigma_c", start_4d, count_4d, mem_dims_4d, norm_sigma_c);

      const int f2 = read_u16_4d_dataset_slice(cont_group,  //
                                               "c_therm_opacity_k_c", start_4d, count_4d, mem_dims_4d, norm_epsilon_c);

      const int f3 = read_u16_4d_dataset_slice(cont_group,  //
                                               "c_therm_emissivity_epsilon_c", start_4d, count_4d, mem_dims_4d,
                                               norm_c_therm_emissivity_epsilon_c);

      if (f1 < 0 || f2 < 0 || f3 < 0) {
        free(norm_sigma_c);
        free(norm_epsilon_c);
        free(norm_c_therm_emissivity_epsilon_c);
        H5Gclose(cont_group);
        return -1;
      }
    }

    for (size_t cell = 0; cell < block_size; cell++) {
      const size_t offset = cell * (size_t)N_frequencies;
      denormalize_continuum_data(cont_to_fill[cell],  //
                                 &norm_sigma_c[offset], &norm_epsilon_c[offset],
                                 &norm_c_therm_emissivity_epsilon_c[offset],                                    //
                                 N_frequencies,                                                                   //
                                 &sigma_c[offset], &epsilon_c[offset], &c_therm_emissivity_epsilon_c[offset]);  //
    }

    free(norm_sigma_c);
    free(norm_epsilon_c);
    free(norm_c_therm_emissivity_epsilon_c);
  }

  H5Gclose(cont_group);
  return 0;
}

/////////////////////////////////////////////////////////////////////////
// continuum_main_demo
/////////////////////////////////////////////////////////////////////////
int
continuum_main_demo(int argc, char **argv) {
  const char *filename = "cont_atmosphere.h5";
  if (argc > 1 && argv[1] != NULL && argv[1][0] != '\0') {
    filename = argv[1];
  }

  hid_t file = H5I_INVALID_HID;
  H5E_BEGIN_TRY { file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT); }
  H5E_END_TRY;
  if (file < 0) {
    if (errno == ENOENT) {
      // Create the file on first run so the continuum demo works out-of-the-box.
      file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (file < 0) {
        fprintf(stderr, "Error creating HDF5 file '%s'\n", filename);
        return EXIT_FAILURE;
      }
      printf("Created new HDF5 file '%s'\n", filename);
    } else {
      fprintf(stderr, "Error opening HDF5 file '%s'\n", filename);
      return EXIT_FAILURE;
    }
  }

  int N_x = 10, N_y = 10, N_z = 5, N_frequencies = 100;

  make_continuum_hdf5_mpi(file, N_x, N_y, N_z, N_frequencies);  // Example dimensions
  THDF_create_continuum_4D_matrices(file, N_x, N_y, N_z, N_frequencies);

  double *sigma_c                        = (double *)malloc(N_frequencies * sizeof(double));
  double *epsilon_c                      = (double *)malloc(N_frequencies * sizeof(double));
  double *c_therm_emissivity_epsilon_c = (double *)malloc(N_frequencies * sizeof(double));

  u_int16_t *norm_sigma_c                        = (uint16_t *)malloc(N_frequencies * sizeof(uint16_t));
  u_int16_t *norm_epsilon_c                      = (uint16_t *)malloc(N_frequencies * sizeof(uint16_t));
  u_int16_t *norm_c_therm_emissivity_epsilon_c = (uint16_t *)malloc(N_frequencies * sizeof(uint16_t));

  for (int i = 0; i < N_frequencies; i++) {
    sigma_c[i]                        = 1.0e-2 + 1.0e-4 * i;
    epsilon_c[i]                      = 1.0e-3 + 1.0e-5 * i;
    c_therm_emissivity_epsilon_c[i] = 1.0e-1 + 1.0e-3 * i;
  }

  THDF_continuum_t cont = normalize_continuum_data(sigma_c, epsilon_c, c_therm_emissivity_epsilon_c, N_frequencies,
                                                   norm_sigma_c, norm_epsilon_c, norm_c_therm_emissivity_epsilon_c);

  write_normalized_continuum_to_hdf5(file, N_frequencies, 0, 0, 0, 1, 1, 1, &cont, norm_sigma_c, norm_epsilon_c,
                                     norm_c_therm_emissivity_epsilon_c);
  write_normalized_continuum_to_hdf5(file, N_frequencies, 1, 1, 1, 1, 1, 1, &cont, norm_sigma_c, norm_epsilon_c,
                                     norm_c_therm_emissivity_epsilon_c);
  write_normalized_continuum_to_hdf5(file, N_frequencies, 2, 2, 2, 1, 1, 1, &cont, norm_sigma_c, norm_epsilon_c,
                                     norm_c_therm_emissivity_epsilon_c);

  H5Fclose(file);

  free(sigma_c);
  free(epsilon_c);
  free(c_therm_emissivity_epsilon_c);

  free(norm_sigma_c);
  free(norm_epsilon_c);
  free(norm_c_therm_emissivity_epsilon_c);

  return EXIT_SUCCESS;
}
