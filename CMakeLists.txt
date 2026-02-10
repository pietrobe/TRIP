cmake_minimum_required(VERSION 3.5...3.19)

set(SOLAR_3D_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${SOLAR_3D_ROOT_PATH}/cmake")

project(par_moonolith_usage VERSION 1.0.0)

### DEBUG OPTION
option(DEBUG "Build with debugging symbols and flags" OFF)
if(DEBUG)
    message(STATUS "Configuring Debug build")
    add_definitions(-DDEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -Wextra -std=c++17")
    # set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -Wall -Wextra -fsanitize=address,undefined -fno-omit-frame-pointer")
    # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address,undefined")
    set(EXECUTABLE_NAME "solar_3D_debug")
else()
    message(STATUS "Configuring Release build")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -std=c++17 -Wno-deprecated-declarations -Wno-tautological-constant-compare -Wno-recommended-option")
    #set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Ofast -fno-signed-zeros -fno-trapping-math -fassociative-math -march=native -mtune=native -mprefer-vector-width=512 -fopenmp -std=c++17 -Wno-deprecated-declarations -Wno-tautological-constant-compare -Wno-recommended-option -pedantic ")
    set(EXECUTABLE_NAME "solar_3D")
endif()

### DEPENDENCIES
find_package(MPI REQUIRED) 
find_package(HDF5 REQUIRED)
find_package(HDF5 COMPONENTS C HL REQUIRED)

if(NOT HDF5_INCLUDE_DIRS)
    set(HDF5_INCLUDE_DIRS "${HDF5_ROOT}/include")
    set(HDF5_LIBRARIES "${HDF5_ROOT}/lib/libhdf5.so;${HDF5_ROOT}/lib/libhdf5_hl.so")
endif()

include(CMakeFindDependencyMacro)
find_package(PkgConfig REQUIRED)

find_package(yaml-cpp QUIET)

if(NOT yaml-cpp_FOUND)
    if (DEFINED YAML_CPP_URL)
        message(STATUS "yaml-cpp not found on system, fetching from URL...")
        include(FetchContent)
        FetchContent_Declare(
            yaml-cpp
            URL ${YAML_CPP_URL}
            # Optional: add hash for integrity check
            # URL_HASH SHA256=<hash_value>
        )

    else()

        message(STATUS "yaml-cpp not found on system, fetching from Git...")
        include(FetchContent)
        FetchContent_Declare(
            yaml-cpp
            GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
            GIT_TAG 0.8.0 # Can be a tag (yaml-cpp-x.x.x), a commit hash, or a branch name (master)
        )
    endif()
    FetchContent_MakeAvailable(yaml-cpp)
else()
    message(STATUS "Found yaml-cpp on system: ${yaml-cpp_VERSION}")
endif()

#target_link_libraries(YOUR_LIBRARY PUBLIC yaml-cpp::yaml-cpp) # The library or executable that require yaml-cpp library

find_package(Boost 1.61 COMPONENTS program_options regex serialization filesystem  atomic  REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
set(Boost_LIBRARY_DIRS "${BOOST_ROOT}/lib")


### SOURCE FILES
list(APPEND SOURCES_FILES 
    src/Faddeeva.cpp 
    src/Formal_solver.cpp 
    src/Legendre_rule.cpp 
    # src/main.cpp 
    src/Rotation_matrix.cpp
    src/RT_problem.cpp
    src/RT_problem_IO.cpp
    src/RT_problem_txt_IO.cpp
    src/RT_solver.cpp
    src/Utilities.cpp
    src/RT_utility.cpp
    src/RT_arbitrary_beam.cpp
)

list(APPEND HEADERS_FILES 
    src/Faddeeva.hpp
    src/Formal_solver.hpp
    src/Legendre_rule.hpp
    src/Rotation_matrix.hpp
    src/RT_problem.hpp
    src/RT_solver.hpp
    src/Utilities.hpp
    src/RT_utility.hpp
    src/RT_arbitrary_beam.hpp
    src/RT_types.hpp
)

### RII LINKAGE
set(RII_INCLUDES "${RII_ROOT_PATH}")
message(${RII_INCLUDES})

option(USE_ACC "Use ACC acceleration" OFF)

if(USE_ACC)
    add_definitions(-DUSE_ACC)
    message(STATUS "ACC acceleration enabled")
    set(ACC_SOLAR_3D _ON_)
else()
    message(STATUS "ACC acceleration disabled")
    set(ACC_SOLAR_3D _OFF_)
endif()

include_directories(${RII_INCLUDES}/src)
include_directories(${RII_INCLUDES}/src/app)
include_directories(${RII_INCLUDES}/src/atmosphere_models)
include_directories(${RII_INCLUDES}/src/bench)
include_directories(${RII_INCLUDES}/src/doppler_shift)
include_directories(${RII_INCLUDES}/src/epsilon_alt)
include_directories(${RII_INCLUDES}/src/factory)
include_directories(${RII_INCLUDES}/src/faddeeva)
include_directories(${RII_INCLUDES}/src/fastgl)
include_directories(${RII_INCLUDES}/src/formal_solver)
include_directories(${RII_INCLUDES}/src/gauss_kronrod_adaptive)
include_directories(${RII_INCLUDES}/src/math)
include_directories(${RII_INCLUDES}/src/md_matrix)
include_directories(${RII_INCLUDES}/src/quadrature_intervals)
include_directories(${RII_INCLUDES}/src/rii_algo)
include_directories(${RII_INCLUDES}/src/rii_algo_emission_coefficient_contrib)
include_directories(${RII_INCLUDES}/src/rii_algo_emission_coefficient_fg)
include_directories(${RII_INCLUDES}/src/rii_base)
include_directories(${RII_INCLUDES}/src/rii_frequencies_grid)
include_directories(${RII_INCLUDES}/src/riii_approx_emission_coeffcient)
include_directories(${RII_INCLUDES}/src/riii_exact_emission_coeffcient)
include_directories(${RII_INCLUDES}/src/rii_include)
include_directories(${RII_INCLUDES}/src/rii_maps)
include_directories(${RII_INCLUDES}/src/rii_RII_KKpQ)
include_directories(${RII_INCLUDES}/src/epsilon_alt)
include_directories(${RII_INCLUDES}/src/sphere_lebedev_rule)
include_directories(${RII_INCLUDES}/src/tk_spline)
include_directories(${RII_INCLUDES}/src/rii_aa_emission_coefficient)
include_directories(${RII_INCLUDES}/src/utils)
include_directories(${RII_INCLUDES}/src/rii_aa_emission_coefficient/)
include_directories(${RII_INCLUDES}/src/acc_includes)
if(USE_ACC)
    include_directories(${RII_INCLUDES}/src_mpi/RII_contrib)
    include_directories(${RII_INCLUDES}/src_acc/src)
endif()

### PETSC DEPENDENCY
message(STATUS "Search Petsc")

set(PETSC_TEST_RUNS TRUE)
set(ENV{PKG_CONFIG_PATH} $ENV{PETSC_DIR}/lib/pkgconfig)

pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

message(STATUS "Check Petsc done")

# Configure header file
set(SOLAR_3D_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(SOLAR_3D_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(SOLAR_3D_VERSION_PATCH ${PROJECT_VERSION_PATCH})

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/configure.h.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/configure.h"
    @ONLY
)

# we use one of our examples here
add_executable(${EXECUTABLE_NAME} ${CMAKE_SOURCE_DIR}/src/main.cpp)

target_include_directories(${EXECUTABLE_NAME} PUBLIC ${HDF5_INCLUDE_DIRS})


add_library(rii SHARED IMPORTED)
set_property(TARGET rii PROPERTY IMPORTED_LOCATION "${RII_ROOT_PATH}/build/librii.so")
message( STATUS  "RII path: ${RII_ROOT_PATH}")

# import all the dependencies and targets # removed sgrid::sgrid
target_link_libraries(${EXECUTABLE_NAME} PUBLIC rii)
target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${Boost_LIBRARIES})
target_link_libraries(${EXECUTABLE_NAME} PUBLIC MPI::MPI_CXX) # Added Melanie
target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${HDF5_LIBRARIES})
target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${HDF5_HL_LIBRARIES})


target_sources(${EXECUTABLE_NAME} PUBLIC ${SOURCES_FILES})
target_include_directories(${EXECUTABLE_NAME} PUBLIC  
    src
    ${CMAKE_CURRENT_BINARY_DIR}  # Add binary dir for configure.h
)


if(NOT yaml-cpp_FOUND)
    target_link_libraries(${EXECUTABLE_NAME} PUBLIC yaml-cpp::yaml-cpp)
else()
    target_link_libraries(${EXECUTABLE_NAME} PUBLIC yaml-cpp)
endif()

target_include_directories(${EXECUTABLE_NAME} SYSTEM PUBLIC ${PETSC_INCLUDES})
target_link_libraries(${EXECUTABLE_NAME} PUBLIC PkgConfig::PETSC)
target_compile_definitions(${EXECUTABLE_NAME} PUBLIC ${PETSC_DEFINITIONS})

# Build thdf/libhdf_output.a using its Makefile
if(TRUE)
    # Collect all C source and header files in src/thdf
    file(GLOB THDF_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/src/thdf/*.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/thdf/*.h"
    )

    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/thdf/libhdf_output.a
        COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR}/src/thdf make CC=mpicc lib
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/thdf
        BYPRODUCTS ${CMAKE_CURRENT_SOURCE_DIR}/src/thdf/libhdf_output.a
        DEPENDS ${THDF_SOURCES}
        COMMENT "Building libhdf_output.a in thdf using Makefile"
    )

    add_custom_target(hdf_output_lib ALL
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/thdf/libhdf_output.a
    )

    # Ensure solar_3D depends on the static library
    add_dependencies(${EXECUTABLE_NAME} hdf_output_lib)

    # Include thdf headers
    target_include_directories(${EXECUTABLE_NAME} PUBLIC src/thdf/)

    # Link the static library
    target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/thdf/libhdf_output.a)

    # Provide a clean target for the thdf subproject
    add_custom_target(thdf_clean
        COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR}/src/thdf make clean
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/thdf
        COMMENT "Cleaning thdf subproject (make clean)"
    )

    # Convenience target to clean both top-level and thdf subproject
    add_custom_target(clean-all
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
        COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR}/src/thdf make clean
        COMMENT "Clean top-level build and thdf subproject"
    )
endif()

# Add Grid Manager
add_subdirectory(src/GridManager)
target_link_libraries(${EXECUTABLE_NAME} PUBLIC GridManager)

