# HDF5 Output Module - CMake Build Guide

This CMakeLists.txt builds the HDF5 output library module for standalone use or as a subproject.

## What Gets Built

- **libhdf_output.a** - Static library (`.a` file) containing all output utilities
- **hdf_output_example** - Example executable demonstrating basic usage
- **hdf_output_mpi_example** - MPI example executable

## Quick Start - Standalone Build

```bash
# From the output directory or its parent
mkdir build
cd build
cmake ..
cmake --build . --config Release
```

This produces:
- `build/libhdf_output.a` - The static library
- `build/hdf_output_example` - Example program
- `build/hdf_output_mpi_example` - MPI example program

## Using as a Subproject

When this module is used as a subproject in another CMake project:

### Method 1: Direct with add_subdirectory()

In your parent `CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.15)
project(my_project)

# Include the hdf_output module as a subproject
add_subdirectory(path/to/hdf5_demo/output)

# Link library to your target
add_executable(my_app src/main.c)
target_link_libraries(my_app PRIVATE hdf_output)
```

### Method 2: Modern approach with FetchContent

```cmake
cmake_minimum_required(VERSION 3.14)
project(my_project)

include(FetchContent)

FetchContent_Declare(
    hdf_output
    GIT_REPOSITORY https://github.com/user/hdf5_demo.git
    GIT_TAG main
    SOURCE_SUBDIR output
)

FetchContent_MakeAvailable(hdf_output)

# Use the library
add_executable(my_app src/main.c)
target_link_libraries(my_app PRIVATE hdf_output)
```

## Key Points

✅ **Static Library Only** - When used as a subproject, only the `.a` library file matters  
✅ **Automatic Dependencies** - HDF5 and system libraries are automatically linked  
✅ **No Manual Configuration** - Include paths and compiler flags are handled automatically  
✅ **Executables Optional** - Example programs are built but don't interfere with subproject use  

## Build Options

```bash
# Release build (default, optimized)
cmake -DCMAKE_BUILD_TYPE=Release ..

# Debug build (with debug symbols)
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Custom installation prefix
cmake -DCMAKE_INSTALL_PREFIX=/opt/hdf5_output ..

# Verbose build output
cmake --build . --verbose
```

## Installation

```bash
# Install to system location (requires sudo for /usr/local)
sudo cmake --install . --prefix /usr/local

# Install to user home directory
cmake --install . --prefix ~/.local/hdf5_output

# Install to custom location
cmake --install . --prefix /opt/mylibs
```

## Example: Using installed library

After installing with `cmake --install . --prefix ~/.local/hdf5_output`:

```cmake
# In your CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(my_app)

# Find the installed package
find_package(hdf_output REQUIRED)

add_executable(my_app src/main.c)
target_link_libraries(my_app PRIVATE hdf5_demo::hdf_output)
```

## Troubleshooting

### HDF5 not found

```bash
# If HDF5 is in a non-standard location:
cmake .. -DHDF5_ROOT=/path/to/hdf5

# Try pkg-config
export PKG_CONFIG_PATH=/path/to/hdf5/lib/pkgconfig:$PKG_CONFIG_PATH
cmake ..
```

### Build fails with compiler errors

Ensure you have:
- GCC 7+ or Clang 5+ (for C18 support)
- libhdf5-dev package installed
- Standard build tools (make, cmake)

On Ubuntu/Debian:
```bash
sudo apt-get install cmake gcc libhdf5-dev
```

On macOS:
```bash
brew install cmake hdf5
```

### Want to build only the library (skip executables)?

```bash
cmake --build . --target hdf_output
```

## CMakeLists.txt Structure

```
Total Components:
├── Library (libhdf_output.a)
│   ├── Source files (*.c)
│   ├── Public headers (*.h)
│   └── Dependencies (HDF5, libm)
│
├── Executables (optional)
│   ├── hdf_output_example
│   └── hdf_output_mpi_example
│
└── Installation & Export targets
    ├── Install rules
    └── CMake package config
```

## File Structure

```
output/
├── CMakeLists.txt           # CMake build configuration (this file)
├── *.c, *.h                 # Source and header files
├── Makefile                 # Legacy Makefile (optional, can coexist)
└── README_CMAKE.md          # This guide
```

## Supported Compiler Flags

The CMake automatically applies:
- `-O3` - Optimize for speed
- `-Wall -Wextra` - Compiler warnings
- `-march=native -mtune=native` - CPU-specific optimization

To override:
```bash
cmake .. -DCMAKE_C_FLAGS="-O2 -Wall"
```

## Integration with package managers

After installation, other CMake projects can find it:

```cmake
# For system installation
find_package(hdf_output REQUIRED)

# For custom prefix
find_package(hdf_output REQUIRED CONFIG PATHS ~/.local/hdf5_output/lib/cmake)
```
