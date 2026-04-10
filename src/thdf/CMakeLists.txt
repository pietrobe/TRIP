cmake_minimum_required(VERSION 3.15)
project(hdf_output 
    VERSION 1.0.0
    DESCRIPTION "HDF5 Output Library and Examples"
    LANGUAGES C)

# Set C standard (use C17 for maximum compatibility with HDF5 and MPI)
set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(_THDF_BUILD_EXAMPLES_DEFAULT ON)
else()
    set(_THDF_BUILD_EXAMPLES_DEFAULT OFF)
endif()
option(THDF_BUILD_EXAMPLES "Build thdf example executables" ${_THDF_BUILD_EXAMPLES_DEFAULT})

# Compiler flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -Wextra -march=native -mtune=native")

# Find HDF5 package
find_package(HDF5 REQUIRED COMPONENTS C HL)

# Find MPI (optional, but HDF5 may be compiled with it)
find_package(MPI QUIET)

# Determine HDF5 target names (HDF5 2.0+ uses different naming)
if(TARGET hdf5::hdf5)
    set(HDF5_TARGET hdf5::hdf5)
    set(HDF5_HL_TARGET hdf5::hdf5_hl)
elseif(TARGET HDF5::C)
    set(HDF5_TARGET HDF5::C)
    set(HDF5_HL_TARGET HDF5::HL)
else()
    set(HDF5_TARGET HDF5::HDF5)
    set(HDF5_HL_TARGET HDF5::HL)
endif()

# ==============================================================================
# STATIC LIBRARY TARGET: libhdf_output.a
# ==============================================================================

add_library(hdf_output STATIC
    output_utils.c
    hdf_output_field.c
    hdf_output_ar_field.c
    hdf_output_JKQ.c
    hdf_output_3D_field.c
    hdf_output_example_3D.c
    hdf_fapl_from_env.c
    hdf_output_domain_decomposition_3D.c
    hdf_output_zslice_3D.c
    hdf_output_single_zslice_3D.c
    hdf_output_example_zslice_3D.c
    hdf_output_example_multi_zslice_3D.c
    hdf_output_example_single_zslice_3D.c
)

# Set include directories for the library
target_include_directories(hdf_output
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

# Link HDF5 libraries to the static library
if(MPI_FOUND)
    target_link_libraries(hdf_output
        PUBLIC
            ${HDF5_TARGET}
            ${HDF5_HL_TARGET}
            MPI::MPI_C
            m
    )
else()
    target_link_libraries(hdf_output
        PUBLIC
            ${HDF5_TARGET}
            ${HDF5_HL_TARGET}
            m
    )
endif()

# Set library properties
set_target_properties(hdf_output PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
)

# ==============================================================================
# EXECUTABLE TARGETS (optional, for standalone builds)
# ==============================================================================

if(THDF_BUILD_EXAMPLES)
    # hdf_output_example executable
    add_executable(hdf_output_example
        hdf_output_example.c
    )

    target_link_libraries(hdf_output_example
        PRIVATE
            hdf_output
    )

    # hdf_output_mpi_example executable
    add_executable(hdf_output_mpi_example
        hdf_output_mpi_example.c
    )

    target_link_libraries(hdf_output_mpi_example
        PRIVATE
            hdf_output
    )
endif()

# ==============================================================================
# INSTALLATION CONFIGURATION
# ==============================================================================

# Install the static library
install(
    TARGETS hdf_output
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

# Install public headers
install(
    FILES
        output_utils.h
        hdf_output_field.h
        hdf_output_ar_field.h
        hdf_output_JKQ.h
        hdf_output_3D_field.h
        hdf_output_example_3D.h
        hdf_fapl_from_env.h
        hdf_output_domain_decomposition_3D.h
        hdf_output_zslice_3D.h
        hdf_output_single_zslice_3D.h
        hdf_output_example_zslice_3D.h
        thdf.h
    DESTINATION include
)

# Install executables (optional)
if(THDF_BUILD_EXAMPLES)
    install(
        TARGETS hdf_output_example hdf_output_mpi_example
        RUNTIME DESTINATION bin
        OPTIONAL
    )
endif()

# ==============================================================================
# CMAKE PACKAGE CONFIGURATION (for use as a subproject)
# ==============================================================================

# Export library target so it can be used by parent projects
export(
    TARGETS hdf_output
    NAMESPACE hdf5_demo::
    FILE "${PROJECT_BINARY_DIR}/hdf_outputTargets.cmake"
)

# Create a simple config file for find_package support
file(WRITE "${PROJECT_BINARY_DIR}/hdf_outputConfig.cmake"
    "include(\"\${CMAKE_CURRENT_LIST_DIR}/hdf_outputTargets.cmake\")\n"
)

# ==============================================================================
# BUILD INFORMATION
# ==============================================================================

message(STATUS "=== HDF5 Output Library ===")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "C Standard: ${CMAKE_C_STANDARD}")
message(STATUS "HDF5 Found: ${HDF5_FOUND}")
message(STATUS "HDF5 Version: ${HDF5_VERSION}")
message(STATUS "")
message(STATUS "Targets:")
message(STATUS "  - hdf_output (STATIC LIBRARY) - Core library")
if(THDF_BUILD_EXAMPLES)
    message(STATUS "  - hdf_output_example (EXECUTABLE) - Basic example")
    message(STATUS "  - hdf_output_mpi_example (EXECUTABLE) - MPI example")
endif()
message(STATUS "")
message(STATUS "To use as a subproject:")
message(STATUS "  add_subdirectory(output)")
message(STATUS "  target_link_libraries(my_app PRIVATE hdf_output)")
message(STATUS "============================")
