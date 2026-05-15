#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "hdf_atmos_cpp.hpp"

// Compile with: g++ -std=c++11 -o test_partial_read test_partial_read.cpp hdf_atmos_functions.c -lhdf5 -lhdf5_hl

int
main(int argc, char *argv[]) {
  std::string filename = "/home/simone/cscs_scratch/PORTA/Model_QuietSun_8x8x128_extract.h5";

  if (argc > 1) {
    filename = argv[1];
  }

  try {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
      throw std::runtime_error("Error opening HDF5 file: " + filename);
    }

    std::cout << "==== Testing read_atmos_3D_partial_data ====" << std::endl;
    std::cout << "Reading column at (i=2, j=3) across all heights\n" << std::endl;

    std::vector<THDF_atmos_t> atmos_data;

    int start_i = 2, start_j = 3, start_k = 0;
    int count_i = 1, count_j = 1, count_k = 143;

    int status = hdf_atmos_cpp::read_atmos_3D_partial_data(file, atmos_data, start_i, start_j, start_k, count_i, count_j, count_k);

    if (status < 0) {
      H5Fclose(file);
      throw std::runtime_error("Error reading partial atmosphere data");
    }

    std::cout << "Successfully read " << atmos_data.size() << " data points\n" << std::endl;
    std::cout << std::left << std::setw(8) << "Height(k)"
              << std::setw(20) << "Temperature (K)"
              << std::setw(15) << "Vmicro (km/s)"
              << std::setw(15) << "Damping"
              << std::endl;
    std::cout << std::string(58, '-') << std::endl;

    for (size_t k = 0; k < atmos_data.size(); ++k) {
      const THDF_atmos_t &point = atmos_data[k];
      std::cout << std::left
                << std::setw(8) << k
                << std::setw(20) << std::scientific << std::setprecision(6) << point.temperature
                << std::setw(15) << point.vmicro
                << std::setw(15) << point.damping
                << std::endl;
    }

    H5Fclose(file);
    std::cout << "\nTest completed successfully!" << std::endl;

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
