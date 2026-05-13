#include <iostream>
#include <string>
#include <vector>
#include "hdf_atmos_cpp.hpp"

// Compile with: g++ -std=c++11 -o hdf_atmos_cpp_example hdf_atmos_cpp_example.cpp -lhdf5 -lhdf5_hl

#ifndef HDF_ATMOS_BUILD_LIBRARY
int
main(int argc, char *argv[]) {
  std::string filename = "/home/simone/cscs_scratch/PORTA/Model_QuietSun_8x8x128_extract.h5";

  if (argc > 1) {
    filename = argv[1];
  }

  try {
    // Using the all-in-one convenience function with validation
    hdf_atmos_cpp::validate_and_read_file(filename);

    std::cout << "\n=== Alternative: Using individual read functions ===\n" << std::endl;

    // Alternative: Read data piece by piece
    int                       N_x, N_y, N_z, N_frequencies;
    std::vector<double>       heights;
    std::vector<double>       frequencies;
    THDF_atom_two_levels_t    atom;
    std::vector<THDF_atmos_t> atmos_data;

    hdf_atmos_cpp::read_file(filename, N_x, N_y, N_z, heights, N_frequencies, frequencies, atom, atmos_data);

    std::cout << "Data loaded successfully!" << std::endl;
    std::cout << "Grid dimensions: " << N_x << " x " << N_y << " x " << N_z << std::endl;
    std::cout << "Total atmosphere points: " << atmos_data.size() << std::endl;
    std::cout << "Heights vector size: " << heights.size() << std::endl;
    std::cout << "Frequencies vector size: " << frequencies.size() << std::endl;

    // Access data using std::vector
    std::cout << "\nFirst 5 heights (using std::vector):" << std::endl;
    for (size_t i = 0; i < std::min(size_t(5), heights.size()); ++i) {
      std::cout << "  heights[" << i << "] = " << heights[i] << " cm" << std::endl;
    }

    if (!frequencies.empty()) {
      std::cout << "\nFirst 5 frequencies (using std::vector):" << std::endl;
      for (size_t i = 0; i < std::min(size_t(5), frequencies.size()); ++i) {
        std::cout << "  frequencies[" << i << "] = " << frequencies[i] << " Hz" << std::endl;
      }
    }

    // Access atmosphere data
    std::cout << "\nAtmosphere data at point (5, 10, 3):" << std::endl;
    int idx = 5 * N_y * N_z + 10 * N_z + 3;
    if (idx < static_cast<int>(atmos_data.size())) {
      const THDF_atmos_t &pt = atmos_data[idx];
      std::cout << "  Index: (" << pt.index_i << ", " << pt.index_j << ", " << pt.index_k << ")" << std::endl;
      std::cout << "  Temperature: " << pt.temperature << " K" << std::endl;
      std::cout << "  Vmicro: " << pt.vmicro << " km/s" << std::endl;
    }

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
#endif /* HDF_ATMOS_BUILD_LIBRARY */
