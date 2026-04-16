# HDF5 Atmosphere C++ API

This document describes the C++ inline read functions for HDF5 atmosphere data using `std::vector` for dynamic storage.

## Header File

`hdf_atmos_cpp.h` - Contains all inline C++ read functions in the `hdf_atmos_cpp` namespace.

## Features

- **Inline functions**: All functions are defined in the header file for optimal performance
- **std::vector**: Automatic memory management using C++ standard library
- **Exception handling**: Uses C++ exceptions for error reporting
- **Type safety**: C++ references and const-correctness

## API Functions

### Individual Read Functions

#### `read_geometry()`
```cpp
inline int read_geometry(hid_t file, int &N_x, int &N_y, int &N_z, 
                        std::vector<double> &heights);
```
Reads grid dimensions and height scale from `/geometry` group.

#### `read_frequency_grid()`
```cpp
inline int read_frequency_grid(hid_t file, int &N_frequencies, 
                               std::vector<double> &frequencies);
```
Reads frequency grid data from `/frequency_grid` group.

#### `read_atom()`
```cpp
inline int read_atom(hid_t file, hdf_atom &atom);
```
Reads atomic transition parameters from `/atom` group.

#### `read_atmos_3D_data()`
```cpp
inline int read_atmos_3D_data(hid_t file, std::vector<hdf_atmos> &atmos_data, 
                              int &N_x, int &N_y, int &N_z);
```
Reads 3D atmosphere data from `/atmosphere/atmos_data` dataset. Automatically resizes the vector.

### Convenience Functions

#### `read_file()`
```cpp
inline int read_file(const std::string &filename, 
                    int &N_x, int &N_y, int &N_z,
                    std::vector<double> &heights,
                    int &N_frequencies,
                    std::vector<double> &frequencies,
                    hdf_atom &atom,
                    std::vector<hdf_atmos> &atmos_data);
```
All-in-one function that reads an entire HDF5 file. Throws exceptions on error.

#### `validate_and_read_file()`
```cpp
inline void validate_and_read_file(const std::string &filename);
```
Reads and validates an HDF5 file, printing detailed information about the contents.

## Usage Example

```cpp
#include "hdf_atmos_cpp.h"

int main() {
    // Method 1: Simple validation
    hdf_atmos_cpp::validate_and_read_file("atmosphere.h5");
    
    // Method 2: Read data for processing
    int N_x, N_y, N_z, N_frequencies;
    std::vector<double> heights;
    std::vector<double> frequencies;
    hdf_atom atom;
    std::vector<hdf_atmos> atmos_data;
    
    hdf_atmos_cpp::read_file("atmosphere.h5", N_x, N_y, N_z, heights, 
                             N_frequencies, frequencies, atom, atmos_data);
    
    // Access data using std::vector
    for (size_t i = 0; i < heights.size(); ++i) {
        std::cout << "Height[" << i << "] = " << heights[i] << std::endl;
    }
    
    // Access 3D atmosphere data
    int idx = i * N_y * N_z + j * N_z + k;
    const hdf_atmos &point = atmos_data[idx];
    std::cout << "Temperature: " << point.temperature << std::endl;
    
    return 0;
}
```

## Building

Compile with C++11 or later:

```bash
# Build the C++ example
make cpp

# Run the C++ example
make run-cpp

# Or manually:
g++ -std=c++11 -o program program.cpp -lhdf5 -lhdf5_hl
```

## Advantages over C API

1. **Memory Management**: `std::vector` handles allocation/deallocation automatically
2. **Type Safety**: C++ references prevent null pointer errors
3. **Exception Safety**: RAII ensures resources are properly cleaned up
4. **Convenience**: String handling with `std::string` instead of `char*`
5. **Performance**: Inline functions can be optimized by the compiler
6. **Modern C++**: Compatible with C++11 and later standards

## Compatibility

- Works with files created by the C version (`hdf_atmos_example`)
- Works with files created by the MPI version (`hdf_atmos_mpi_example`)
- Compatible with HDF5 1.14.x library
- Requires C++11 or later compiler

## Error Handling

Functions return:
- `0` on success
- `-1` on error (for individual functions)
- Throw `std::runtime_error` exceptions (for `read_file()` and `validate_and_read_file()`)

## Memory Usage

`std::vector` automatically manages memory:
- `heights`: N_z × sizeof(double) bytes
- `frequencies`: N_frequencies × sizeof(double) bytes
- `atmos_data`: N_x × N_y × N_z × sizeof(hdf_atmos) bytes

For a 30×50×20 grid:
- ~30,000 atmosphere points
- ~2.6 MB for atmosphere data (assuming ~90 bytes per point)
