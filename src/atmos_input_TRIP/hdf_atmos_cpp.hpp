#ifndef HDF_ATMOS_CPP_H
#define HDF_ATMOS_CPP_H

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdexcept>
#include <string>
#include <vector>

#include "hdf_atmos.h"

namespace hdf_atmos_cpp {

  //==========================================================================
  // is_valid_hdf5_file
  //==========================================================================
  inline bool
  is_valid_hdf5_file(const std::string &file_name) {
    htri_t result = H5Fis_hdf5(file_name.c_str());
    return result > 0;
  }

  //==========================================================================
  // create_atmos_compound_type
  //==========================================================================
  // Helper function to create HDF5 compound datatype for THDF_atmos_t structure
  inline hid_t
  create_atmos_compound_type() {
    return create_hdf_atmos_compound_type();
  }

  //==========================================================================
  // read_atmos_3D_partial_data
  //==========================================================================
  // Read partial 3D atmosphere data from HDF5 file
  inline int
  read_atmos_3D_partial_data(hid_t file, std::vector<THDF_atmos_t> &atmos_data, int start_i, int start_j, int start_k,
                             int count_i, int count_j, int count_k) {
    hid_t dataset = H5Dopen2(file, "/atmos", H5P_DEFAULT);
    if (dataset < 0) {
      return -1;
    }

    // Get the dataspace
    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0) {
      H5Dclose(dataset);
      return -1;
    }

    // Define hyperslab for partial read
    hsize_t start[3] = {static_cast<hsize_t>(start_i), static_cast<hsize_t>(start_j), static_cast<hsize_t>(start_k)};
    hsize_t count[3] = {static_cast<hsize_t>(count_i), static_cast<hsize_t>(count_j), static_cast<hsize_t>(count_k)};

    if (H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
      H5Sclose(dataspace);
      H5Dclose(dataset);
      return -1;
    }

    // Resize vector to hold partial data
    size_t total_size = count_i * count_j * count_k;
    atmos_data.resize(total_size);

    // Create memory dataspace
    hid_t memspace = H5Screate_simple(3, count, NULL);
    if (memspace < 0) {
      H5Sclose(dataspace);
      H5Dclose(dataset);
      return -1;
    }

    // Create compound datatype
    hid_t atmos_type = create_atmos_compound_type();

    // Read the partial data
    herr_t status = H5Dread(dataset, atmos_type, memspace, dataspace, H5P_DEFAULT, atmos_data.data());

    H5Tclose(atmos_type);
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    return (status < 0) ? -1 : 0;
  }

  //==========================================================================
  // read_geometry
  //==========================================================================
  // Read geometry data from HDF5 file
  inline int
  read_geometry(hid_t file, int &N_x, int &N_y, int &N_z, double &delta, std::vector<double> &heights) {
    hid_t group = H5Gopen2(file, "/geometry", H5P_DEFAULT);
    if (group < 0) {
      return -1;
    }

    if (H5LTread_dataset_int(group, "N_ij", &N_x) < 0 || H5LTread_dataset_int(group, "N_height", &N_z) < 0) {
      H5Gclose(group);
      return -1;
    }
    N_y = N_x;

    if (H5LTread_dataset_double(group, "delta", &delta) < 0) {
      H5Gclose(group);
      return -1;
    }

    heights.resize(N_z);
    if (H5LTread_dataset_double(group, "height", heights.data()) < 0) {
      H5Gclose(group);
      return -1;
    }

    H5Gclose(group);
    return 0;
  }

  //==========================================================================
  // read_frequency_grid
  //==========================================================================
  // Read frequency grid from HDF5 file
  inline int
  read_frequency_grid(hid_t file, int &N_frequencies, std::vector<double> &frequencies) {
    hid_t group = H5Gopen2(file, "/frequency_grid", H5P_DEFAULT);
    if (group < 0) {
      return -1;
    }

    if (H5LTread_dataset_int(group, "N_frequencies", &N_frequencies) < 0) {
      H5Gclose(group);
      return -1;
    }

    frequencies.clear();

    H5Gclose(group);
    return 0;
  }

  //==========================================================================
  // read_atom
  //==========================================================================
  // Read atom data from HDF5 file
  inline int
  read_atom(hid_t file, THDF_atom_two_levels_t &atom) {
    hid_t dataset = H5Dopen2(file, "/atom_data", H5P_DEFAULT);
    if (dataset < 0) {
      return -1;
    }
    hid_t  atom_type = create_atom_compound_type();
    herr_t status    = H5Dread(dataset, atom_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &atom);
    H5Tclose(atom_type);
    H5Dclose(dataset);
    return (status < 0) ? -1 : 0;
  }

  /**
   * @brief Check if the atom data in the HDF5 file is a two-term atom
   */
  inline bool
  is_two_term_atom(hid_t file) {
    return THDF_is_two_term_atom(file);
  }

  /**
   * @brief Read two-term atom data from HDF5 file
   */
  inline int
  read_atom_two_terms(hid_t file, THDF_atom_two_terms_t &atom) {
    return THDF_read_atom_two_terms_from_hdf5(file, &atom);
  }

  //==========================================================================
  // read_atmos_3D_data
  //==========================================================================
  // Read 3D atmosphere data from HDF5 file
  inline int
  read_atmos_3D_data(hid_t file, std::vector<THDF_atmos_t> &atmos_data, int &N_x, int &N_y, int &N_z) {
    hid_t dataset = H5Dopen2(file, "/atmos", H5P_DEFAULT);
    if (dataset < 0) {
      return -1;
    }

    // Get the dataspace to determine dimensions
    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0) {
      H5Dclose(dataset);
      return -1;
    }

    // Read dimensions from dataspace
    hsize_t dims[3];
    if (H5Sget_simple_extent_dims(dataspace, dims, NULL) < 0) {
      H5Sclose(dataspace);
      H5Dclose(dataset);
      return -1;
    }

    N_x = dims[0];
    N_y = dims[1];
    N_z = dims[2];

    // Resize vector to hold all data
    size_t total_size = N_x * N_y * N_z;
    atmos_data.resize(total_size);

    // Create compound datatype for reading
    hid_t atmos_type = create_atmos_compound_type();

    // Read the data
    herr_t status = H5Dread(dataset, atmos_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, atmos_data.data());

    H5Tclose(atmos_type);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    return (status < 0) ? -1 : 0;
  }

  //==========================================================================
  // read_file
  //==========================================================================
  // Read entire file and return all data
  inline int
  read_file(const std::string &filename, int &N_x, int &N_y, int &N_z, std::vector<double> &heights, int &N_frequencies,
            std::vector<double> &frequencies, THDF_atom_two_levels_t &atom, std::vector<THDF_atmos_t> &atmos_data) {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
      throw std::runtime_error("Error opening HDF5 file: " + filename);
    }

    double delta;  // Not used in this function, but read from file for completeness

    try {
      if (read_geometry(file, N_x, N_y, N_z, delta, heights) < 0) {
        H5Fclose(file);
        throw std::runtime_error("Error reading geometry data");
      }

      if (read_frequency_grid(file, N_frequencies, frequencies) < 0) {
        H5Fclose(file);
        throw std::runtime_error("Error reading frequency grid data");
      }

      if (read_atom(file, atom) < 0) {
        H5Fclose(file);
        throw std::runtime_error("Error reading atom data");
      }

      if (read_atmos_3D_data(file, atmos_data, N_x, N_y, N_z) < 0) {
        H5Fclose(file);
        throw std::runtime_error("Error reading atmosphere data");
      }

      H5Fclose(file);
      return 0;
    } catch (...) {
      H5Fclose(file);
      throw;
    }
  }

  //==========================================================================
  // validate_and_read_file
  //==========================================================================
  // Convenience function with validation
  inline void
  validate_and_read_file(const std::string &filename) {
    int                       N_x, N_y, N_z, N_frequencies;
    std::vector<double>       heights;
    std::vector<double>       frequencies;
    THDF_atom_two_levels_t    atom;
    std::vector<THDF_atmos_t> atmos_data;

    printf("==== Reading HDF5 file: %s ====\n\n", filename.c_str());

    read_file(filename, N_x, N_y, N_z, heights, N_frequencies, frequencies, atom, atmos_data);

    printf("Geometry: N_x=%d, N_y=%d, N_z=%d\n", N_x, N_y, N_z);
    printf("Heights:\n");
    for (int k = 0; k < N_z; k++) {
      printf("  height[%d] = %e cm\n", k, heights[k]);
    }

    printf("Frequency grid: N_frequencies=%d\n", N_frequencies);
    if (!frequencies.empty()) {
      printf("Frequencies:\n");
      for (int i = 0; i < N_frequencies; i++) {
        printf("  frequency[%d] = %e Hz\n", i, frequencies[i]);
      }
    } else {
      printf("(No frequency data in file)\n");
    }

    printf("Atom data:\n");
    printf("  atomic_mass = %f\n", atom.atomic_mass);
    printf("  E_lower = %f\n", atom.E_lower);
    printf("  E_upper = %f\n", atom.E_upper);
    printf("  g_lower = %f\n", atom.g_lower);
    printf("  g_upper = %f\n", atom.g_upper);
    printf("  jl2 = %d\n", atom.jl2);
    printf("  ju2 = %d\n", atom.ju2);
    printf("  Aul = %e\n", atom.Aul);

    // Validate indices
    for (int k = 0; k < N_z; k++) {
      for (int j = 0; j < N_y; j++) {
        for (int i = 0; i < N_x; i++) {
          int                 idx   = i * N_y * N_z + j * N_z + k;
          const THDF_atmos_t &point = atmos_data[idx];
          if (point.index_i != i || point.index_j != j || point.index_k != k) {
            throw std::runtime_error("Data index mismatch at (" + std::to_string(i) + ", " + std::to_string(j) + ", " +
                                     std::to_string(k) + ")");
          }
        }
      }
    }

    printf("Atmosphere data point at (i=0, j=0, k=0):\n");
    const THDF_atmos_t &point = atmos_data[33];
    printf("  temperature = %f\n", point.temperature);
    printf("  vmicro = %f\n", point.vmicro);
    printf("  Doppler_width = %f\n", point.Doppler_width);
    printf("  damping = %f\n", point.damping);
    printf("  pop_lower_level = %f\n", point.pop_lower_level);
    printf("  pop_upper_level = %f\n", point.pop_upper_level);
    printf("  Cul = %f\n", point.Cul);
    printf("  Qel = %f\n", point.Qel);
    printf("  nH = %f\n", point.nH);
    printf("  magnetic_field_x = %f\n", point.magnetic_field_x);
    printf("  magnetic_field_y = %f\n", point.magnetic_field_y);
    printf("  magnetic_field_z = %f\n", point.magnetic_field_z);
    printf("  bulk_velocity_x = %f\n", point.bulk_velocity_x);
    printf("  bulk_velocity_y = %f\n", point.bulk_velocity_y);
    printf("  bulk_velocity_z = %f\n", point.bulk_velocity_z);
    printf("  c_scat_opacity_sigma_c = %f\n", point.c_scat_opacity_sigma_c);
    printf("  c_therm_opacity_k_c = %f\n", point.c_therm_opacity_k_c);
    printf("  c_therm_emissivity_epsilon_cth = %f\n", point.c_therm_emissivity_epsilon_c);

    printf("\nFile validation successful!\n");
  }

}  // namespace hdf_atmos_cpp

#endif  // HDF_ATMOS_CPP_H
