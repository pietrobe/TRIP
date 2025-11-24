cmake_minimum_required(VERSION 3.5...3.19)

set(SOLAR_3D_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${SOLAR_3D_ROOT_PATH}/cmake")

project(par_moonolith_usage VERSION 1.0.0)

include(CMakeFindDependencyMacro)
find_package(PkgConfig REQUIRED)

find_package(sgrid 0 REQUIRED)
find_package(Boost 1.61 COMPONENTS program_options regex serialization filesystem system atomic  REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
set(Boost_LIBRARY_DIRS "${BOOST_ROOT}/lib")

list(APPEND SOURCES_FILES 
    src/Faddeeva.cpp 
    src/Formal_solver.cpp 
    src/Legendre_rule.cpp 
    src/main.cpp 
    src/Rotation_matrix.cpp
    src/RT_problem.cpp
    src/RT_solver.cpp
    src/Utilities.cpp
    src/Test_rii_include.cpp
    src/RT_utility.cpp
)

list(APPEND HEADERS_FILES 
    src/Faddeeva.hpp
    src/Formal_solver.hpp
    src/Legendre_rule.hpp
    src/Rotation_matrix.hpp
    src/RT_problem.hpp
    src/RT_solver.hpp
    src/Utilities.hpp
    src/Test_rii_include.hpp
    src/RT_utility.hpp
)

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

# find dependencies
message(STATUS "Search Petsc")

set(PETSC_TEST_RUNS TRUE)
set(ENV{PKG_CONFIG_PATH} $ENV{PETSC_DIR}/lib/pkgconfig)

pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

message(STATUS "Check Petsc done")

# if(PETSC_FOUND)
#     message(STATUS "Enter Petsc")
#     list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${PETSC_INCLUDES})
#     list(APPEND UTOPIA_THIRDPARTY_LIBRARIES ${PETSC_LIBRARIES})

#     set(CMAKE_C_COMPILER ${PETSC_COMPILER})

#     if(NOT MPI_CXX_COMPILER)
#         set(MPI_CXX_COMPILER $ENV{MPI_CXX_COMPILER})
#         message(STATUS "compiler ${MPI_CXX_COMPILER}")
#     endif()

#     if(MPI_CXX_COMPILER)
#         set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
#         set(CMAKE_CXX_COMPILER_DEBUG ${MPI_CXX_COMPILER})
#     else()
#         execute_process(COMMAND mpicxx -v RESULT_VARIABLE MPICXX_FAILED)

#         if(MPICXX_FAILED)
#             message(
#                 STATUS
#                     "Using CMAKE compiler, you can define MPI_CXX_COMPILER=<alias_or_path_to_your_compiler>"
#             )
#         else()
#             message(STATUS "-----------------------------------------------")
#             message(
#                 STATUS
#                     "\n[MPI] using mpicxx for compiling c++ files.\nIf you want to use your own compiler define MPI_CXX_COMPILER=<alias_or_path_to_your_compiler>"
#             )
#             message(STATUS "-----------------------------------------------")
#             set(CMAKE_CXX_COMPILER mpicxx)
#             set(CMAKE_CXX_COMPILER_DEBUG mpicxx)
#         endif()
#     endif()

#     set(UTOPIA_WITH_PETSC TRUE)
# else()
#     message(WARNING "[Warning] Petsc not found")
#     set(UTOPIA_WITH_PETSC FALSE)
# endif()

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
add_executable(solar_3D ${CMAKE_SOURCE_DIR}/src/main.cpp)
set(CMAKE_CXX_FLAGS  " -O3 -fopenmp -std=c++17 -Wno-deprecated-declarations -Wno-tautological-constant-compare -Wno-recommended-option ")

#set(CMAKE_CXX_FLAGS  " -Ofast -fno-signed-zeros -fno-trapping-math -fassociative-math -march=native -mtune=native -mprefer-vector-width=512 -fopenmp -std=c++17 -Wno-deprecated-declarations -Wno-tautological-constant-compare -Wno-recommended-option -pedantic ")

add_library(rii SHARED IMPORTED)
set_property(TARGET rii PROPERTY IMPORTED_LOCATION "${RII_ROOT_PATH}/build/librii.so")

message( STATUS  "RII path: ${RII_ROOT_PATH}")

# import all the dependencies and targets
target_link_libraries(solar_3D PUBLIC  sgrid::sgrid rii ${Boost_LIBRARIES})

target_sources(solar_3D PUBLIC ${SOURCES_FILES})
target_include_directories(solar_3D PUBLIC  
    src
    ${CMAKE_CURRENT_BINARY_DIR}  # Add binary dir for configure.h
)

target_include_directories(solar_3D SYSTEM PUBLIC ${PETSC_INCLUDES})
# target_link_libraries(solar_3D PUBLIC ${PETSC_LIBRARIES})
target_link_libraries(solar_3D PUBLIC PkgConfig::PETSC)
target_compile_definitions(solar_3D PUBLIC ${PETSC_DEFINITIONS})


