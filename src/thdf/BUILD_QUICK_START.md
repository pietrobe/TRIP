# Quick Build Reference - HDF5 Output Module

## Standard Build (Standalone)

```bash
cd /home/simone/git/hdf5_demo/output
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j$(nproc)
```

**Result:**
- `build/libhdf_output.a` - Static library
- `build/hdf_output_example` - Example program  
- `build/hdf_output_mpi_example` - MPI example

## Using as Subproject in CMake

Your parent project's `CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.15)
project(my_project)

# Add hdf_output as subproject
add_subdirectory(hdf5_demo/output)

# Link your app
add_executable(my_app main.c)
target_link_libraries(my_app PRIVATE hdf_output)
```

**That's it!** Everything (headers, HDF5 libs) is automatic.

## Installation

```bash
# System installation
sudo cmake --install build --prefix /usr/local

# User home
cmake --install build --prefix ~/.local/hdf5_output
```

## What Gets Built

| Target | Type | Description |
|--------|------|-------------|
| **hdf_output** | Static Library (.a) | Core output utilities |
| hdf_output_example | Executable | Basic demo |
| hdf_output_mpi_example | Executable | MPI demo |

## Key Files

- `CMakeLists.txt` - Build configuration
- `README_CMAKE.md` - Full documentation
- `CMakeLists_parent_example.txt` - Parent project template
- `example_subproject_usage.c` - Sample code

## HDF5 Dependency

Ensure HDF5 is installed:

**Ubuntu/Debian:**
```bash
sudo apt-get install libhdf5-dev
```

**macOS:**
```bash
brew install hdf5
```

**CentOS/RHEL:**
```bash
sudo yum install hdf5-devel
```

## Troubleshooting

**CMake can't find HDF5:**
```bash
cmake .. -DHDF5_ROOT=/path/to/hdf5
```

**Want library only (skip executables)?**
```bash
cmake --build build --target hdf_output
```

**Clean rebuild:**
```bash
rm -rf build && mkdir build && cd build && cmake ..
```
