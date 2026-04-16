cmake_minimum_required(VERSION 3.15)
project(hdf_atmos C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)

# Find HDF5 (C + HL components)
find_package(HDF5 REQUIRED COMPONENTS C HL)

# ----------------------------------------------------------------------------
# Static library: hdf_atmos
# hdf_atmos_cpp_example.cpp is intentionally excluded (example only).
# ----------------------------------------------------------------------------
set(LIB_SOURCES
    hdf_atmos.c
    hdf_atmos_functions.c
    hdf_atmos_continuum.c
    hdf_atmos_example_trip.c
)

add_library(hdf_atmos STATIC ${LIB_SOURCES})

# Define the macro that suppresses main() symbols in guarded source files
target_compile_definitions(hdf_atmos PRIVATE HDF_ATMOS_BUILD_LIBRARY)

target_include_directories(hdf_atmos
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${HDF5_INCLUDE_DIRS}
)

target_link_libraries(hdf_atmos
    PUBLIC
        ${HDF5_C_LIBRARIES}
        ${HDF5_HL_LIBRARIES}
        m
)

# ----------------------------------------------------------------------------
# Optional MPI support
# cmake -DHDF_ATMOS_WITH_MPI=ON ...
# ----------------------------------------------------------------------------
option(HDF_ATMOS_WITH_MPI "Build with MPI support (adds hdf_atmos_mpi.c)" OFF)

if(HDF_ATMOS_WITH_MPI)
    find_package(MPI REQUIRED COMPONENTS C)
    target_sources(hdf_atmos PRIVATE hdf_atmos_mpi.c)
    target_link_libraries(hdf_atmos PUBLIC MPI::MPI_C)
endif()

# ----------------------------------------------------------------------------
# Install rules
# ----------------------------------------------------------------------------
install(TARGETS hdf_atmos
    ARCHIVE DESTINATION lib
)

install(FILES
    hdf_atmos.h
    hdf_atmos_continuum.h
    hdf_atmos_cpp.hpp
    DESTINATION include
)
