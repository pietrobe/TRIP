# Building and Using libhdf_atmos

This document covers how to build `libhdf_atmos.a` with CMake and how to link it from an external project.

## Requirements

- CMake 3.15+
- C99 / C++11 compiler
- HDF5 library with the HL component (`libhdf5`, `libhdf5_hl`)
- *(optional)* MPI for the parallel I/O sources

## Build the library

```bash
# from the repo root or directly from this directory
cmake -S atmos_input_TRIP -B atmos_input_TRIP/build
cmake --build atmos_input_TRIP/build
```

The archive is written to `atmos_input_TRIP/build/libhdf_atmos.a`.

### Options

| Option | Default | Description |
|---|---|---|
| `HDF_ATMOS_WITH_MPI` | `OFF` | Include `hdf_atmos_mpi.c` and link against MPI |

```bash
cmake -S atmos_input_TRIP -B atmos_input_TRIP/build -DHDF_ATMOS_WITH_MPI=ON
```

### Install to a prefix

```bash
cmake -S atmos_input_TRIP -B atmos_input_TRIP/build \
      -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --install atmos_input_TRIP/build
```

This installs:
```
/path/to/install/
  lib/
    libhdf_atmos.a
  include/
    hdf_atmos.h
    hdf_atmos_continuum.h
    hdf_atmos_cpp.hpp
```

---

## Linking from an external CMake project

### Method 1 — `add_subdirectory`

If you keep `atmos_input_TRIP` inside your project tree:

```cmake
# CMakeLists.txt of the consuming project
add_subdirectory(path/to/atmos_input_TRIP)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE hdf_atmos)
```

`target_link_libraries` automatically propagates the include directories and
the HDF5 link flags — no extra `target_include_directories` call is needed.

### Method 2 — installed prefix (`find_package` / manual)

After running `cmake --install`, point your project to the prefix:

```cmake
# Locate headers and archive manually
set(HDF_ATMOS_PREFIX /path/to/install)

add_executable(my_app main.cpp)

target_include_directories(my_app PRIVATE ${HDF_ATMOS_PREFIX}/include)
target_link_libraries(my_app PRIVATE
    ${HDF_ATMOS_PREFIX}/lib/libhdf_atmos.a
    hdf5 hdf5_hl m
)
```

---

## Including headers

For the **C API** (pure C or mixed C/C++ projects):

```c
#include "hdf_atmos.h"
#include "hdf_atmos_continuum.h"
```

For the **C++ wrapper** (header-only, `hdf_atmos_cpp` namespace):

```cpp
#include "hdf_atmos_cpp.hpp"
```

The C++ header wraps all functions inline and requires no additional compilation unit.

---

## `HDF_ATMOS_BUILD_LIBRARY` macro

The source files `hdf_atmos_mpi.c` and `hdf_atmos_cpp_example.cpp` each contain
a `main()` function gated behind this macro:

```c
#ifndef HDF_ATMOS_BUILD_LIBRARY
int main(...) { ... }
#endif
```

CMake defines `HDF_ATMOS_BUILD_LIBRARY` automatically when compiling the
library target, so no `main` symbol is ever emitted into the archive.
If you compile those files by hand and want to suppress `main`, pass
`-DHDF_ATMOS_BUILD_LIBRARY` to the compiler.

---

## Minimal usage example (C++)

```cpp
#include "hdf_atmos_cpp.hpp"
#include <iostream>
#include <vector>

int main() {
    const std::string file = "atmosphere_mpi.h5";

    // All-in-one: reads and prints the file contents
    hdf_atmos_cpp::validate_and_read_file(file);

    // Or read into your own variables
    int N_x, N_y, N_z, N_frequencies;
    std::vector<double>       heights, frequencies;
    THDF_atom_two_levels_t    atom;
    std::vector<THDF_atmos_t> atmos_data;

    hdf_atmos_cpp::read_file(file, N_x, N_y, N_z,
                             heights, N_frequencies, frequencies,
                             atom, atmos_data);

    std::cout << "Grid: " << N_x << " x " << N_y << " x " << N_z << "\n";

    // Index a point: (i, j, k)
    int i = 2, j = 3, k = 1;
    const THDF_atmos_t &pt = atmos_data[i * N_y * N_z + j * N_z + k];
    std::cout << "T = " << pt.temperature << " K\n";

    return 0;
}
```

Compile manually (without CMake):

```bash
g++ -std=c++11 -o my_app main.cpp \
    -I/path/to/install/include \
    /path/to/install/lib/libhdf_atmos.a \
    -lhdf5 -lhdf5_hl -lm
```
