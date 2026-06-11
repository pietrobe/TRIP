cmake_minimum_required(VERSION 3.5...3.19)

set(SOLAR_3D_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${SOLAR_3D_ROOT_PATH}/cmake")

project(par_moonolith_usage VERSION 1.0.0 LANGUAGES C CXX)

# for CTest enabling
# enable_testing()

### DEBUG OPTION
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
	set(OPENMP_FLAG "-qopenmp")
	add_compile_options(-Wno-tautological-constant-compare)
	message(STATUS
		"IntelLLVM compiler detected: using -qopenmp, disabling -Wtautological-constant-compare"
	)
else()
	set(OPENMP_FLAG "-fopenmp")
endif()

# Architecture-specific optimization flags. The ARM and x86 branches are kept
# strictly separate: -mcpu is the canonical target flag on aarch64 (where
# -mtune alone does not enable the native ISA, e.g. NEON/SVE), while on x86 it
# is deprecated/ambiguous and -mtune/-march are used instead. No x86 flag may
# leak into the ARM branch or vice versa.
include(CheckCXXCompilerFlag)
string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" TRIP_SYSTEM_PROCESSOR)
if(TRIP_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)$")
	# ARM: target NVIDIA Grace (Neoverse V2, as in GH200 nodes). Prefer the
	# hardcoded neoverse-v2 over -mcpu=native so that building on a front-end
	# node of a different ARM arch still targets the compute nodes.
	if(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC|PGI")
		# NVHPC does not understand -mcpu; its target-processor flag is -tp.
		set(ARCH_FLAGS "-tp=neoverse-v2")
	else()
		check_cxx_compiler_flag("-mcpu=neoverse-v2" COMPILER_SUPPORTS_NEOVERSE_V2)
		if(COMPILER_SUPPORTS_NEOVERSE_V2)
			set(ARCH_FLAGS "-mcpu=neoverse-v2")
			# Neoverse V2 implements SVE2 with 128-bit vectors: fixing the
			# scalable vector length lets the compiler emit fixed-width SVE
			# code instead of vector-length-agnostic code.
			set(CMAKE_REQUIRED_FLAGS "-mcpu=neoverse-v2")
			check_cxx_compiler_flag("-msve-vector-bits=128" COMPILER_SUPPORTS_SVE_VECTOR_BITS)
			unset(CMAKE_REQUIRED_FLAGS)
			if(COMPILER_SUPPORTS_SVE_VECTOR_BITS)
				set(ARCH_FLAGS "${ARCH_FLAGS} -msve-vector-bits=128")
			endif()
		else()
			set(ARCH_FLAGS "-mcpu=native")
		endif()
	endif()
	message(STATUS
		"ARM CPU detected (${CMAKE_SYSTEM_PROCESSOR}): targeting Neoverse V2 / Grace"
	)
elseif(TRIP_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|i[3-6]86)$")
	# x86: probe AVX-512 and AVX2 at configure time, then select the best flags.
	# Compiler-flag acceptance is checked first; on Linux the running CPU is also
	# verified via /proc/cpuinfo so we never emit instructions that will SIGILL.
	set(ARCH_FLAGS "")
	set(TRIP_X86_ISA "baseline")

	check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512F)
	if(COMPILER_SUPPORTS_AVX512F)
		set(_cpu_avx512 0)
		if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
			execute_process(
				COMMAND grep -c "avx512f" /proc/cpuinfo
				OUTPUT_VARIABLE _cpu_avx512
				OUTPUT_STRIP_TRAILING_WHITESPACE
				ERROR_QUIET
			)
		else()
			set(_cpu_avx512 1) # non-Linux: trust the compiler flag check
		endif()
		if(_cpu_avx512 GREATER 0)
			# -march=native auto-enables all AVX-512 sub-extensions present on
			# this CPU (BW, DQ, VL, VNNI, VBMI2, …).
			set(ARCH_FLAGS "-march=native")
			# -mprefer-vector-width=512 tells the auto-vectoriser to use 512-bit
			# registers. Beneficial on modern CPUs (Ice Lake, Sapphire Rapids)
			# that do not throttle clock speed on 512-bit operations.
			check_cxx_compiler_flag("-mprefer-vector-width=512"
				COMPILER_SUPPORTS_PREFER_VECTOR_WIDTH_512
			)
			if(COMPILER_SUPPORTS_PREFER_VECTOR_WIDTH_512)
				set(ARCH_FLAGS "${ARCH_FLAGS} -mprefer-vector-width=512")
			endif()
			set(TRIP_X86_ISA "AVX-512")
		endif()
	endif()

	if(NOT ARCH_FLAGS)
		check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
		set(_cpu_avx2 0)
		if(COMPILER_SUPPORTS_AVX2)
			if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
				execute_process(
					COMMAND grep -c "avx2" /proc/cpuinfo
					OUTPUT_VARIABLE _cpu_avx2
					OUTPUT_STRIP_TRAILING_WHITESPACE
					ERROR_QUIET
				)
			else()
				set(_cpu_avx2 1) # non-Linux: trust the compiler flag check
			endif()
		endif()
		if(_cpu_avx2 GREATER 0)
			set(ARCH_FLAGS "-march=native")
			set(TRIP_X86_ISA "AVX2")
		else()
			set(ARCH_FLAGS "-mtune=native")
			set(TRIP_X86_ISA "baseline (no AVX2/AVX-512)")
		endif()
	endif()

	message(STATUS
		"x86 CPU detected (${CMAKE_SYSTEM_PROCESSOR}): ISA=${TRIP_X86_ISA}"
	)
else()
	set(ARCH_FLAGS "")
	message(WARNING
		"Unknown CPU architecture '${CMAKE_SYSTEM_PROCESSOR}': no architecture-specific flags set"
	)
endif()
message(STATUS "Architecture flags (${CMAKE_SYSTEM_PROCESSOR}): ${ARCH_FLAGS}")

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
	# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -std=c++17 -Wno-deprecated-declarations -Wno-tautological-constant-compare -Wno-recommended-option")
	set(CMAKE_CXX_FLAGS
		"${CMAKE_CXX_FLAGS} -O3 ${OPENMP_FLAG} -std=c++17 -Wno-deprecated-declarations -ffast-math ${ARCH_FLAGS}"
	)
	# set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Ofast -fno-signed-zeros -fno-trapping-math -fassociative-math -march=native -mtune=native -mprefer-vector-width=512 -fopenmp -std=c++17 -Wno-deprecated-declarations -Wno-tautological-constant-compare -Wno-recommended-option -pedantic ")
	set(EXECUTABLE_NAME "solar_3D")
endif()

message(STATUS
	"C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} (${CMAKE_CXX_COMPILER})"
)
message(STATUS "C++ compiler options: ${CMAKE_CXX_FLAGS}")

### DEPENDENCIES
# Skip MPI compiler test which doesn't work with NVHPC compilers
# The test fails but we can link the library directly
set(MPI_SKIP_COMPILER_TEST TRUE)

# Find only the C component. Some environments do not provide a working
# MPI::MPI_CXX target for C++ compilation, but the C target is sufficient for
# linking code that uses the C MPI interface.
find_package(MPI REQUIRED COMPONENTS C)

if(TARGET MPI::MPI_CXX)
	set(TRIP_MPI_TARGET MPI::MPI_CXX)
elseif(TARGET MPI::MPI_C)
	set(TRIP_MPI_TARGET MPI::MPI_C)
	message(STATUS "Using MPI::MPI_C for C++ targets")
else()
	message(FATAL_ERROR "No usable MPI imported target was found")
endif()

find_package(HDF5 REQUIRED)
find_package(HDF5 COMPONENTS C HL REQUIRED)

if(NOT HDF5_INCLUDE_DIRS)
	set(HDF5_INCLUDE_DIRS "${HDF5_ROOT}/include")
	set(HDF5_LIBRARIES
		"${HDF5_ROOT}/lib/libhdf5.so;${HDF5_ROOT}/lib/libhdf5_hl.so"
	)
endif()

include(CMakeFindDependencyMacro)
find_package(PkgConfig REQUIRED)
find_package(yaml-cpp QUIET)

if(NOT yaml-cpp_FOUND)
	if(DEFINED YAML_CPP_URL)
		message(STATUS "yaml-cpp not found on system, fetching from URL...")
		include(FetchContent)
		FetchContent_Declare(
			yaml-cpp
			URL
			${YAML_CPP_URL}
			# Optional: add hash for integrity check
			# URL_HASH SHA256=<hash_value>
		)

	else()

		message(STATUS "yaml-cpp not found on system, fetching from Git...")
		include(FetchContent)
		FetchContent_Declare(
			yaml-cpp
			GIT_REPOSITORY
			https://github.com/jbeder/yaml-cpp.git
			GIT_TAG
			0.8.0
			# Can be a tag (yaml-cpp-x.x.x), a commit hash, or a branch name (master)
		)
	endif()
	FetchContent_MakeAvailable(yaml-cpp)
else()
	message(STATUS "Found yaml-cpp on system: ${yaml-cpp_VERSION}")
endif()

find_package(
	Boost
	1.61
	COMPONENTS
		program_options
		regex
		serialization
		filesystem
		atomic
	REQUIRED
)
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
	src/RT_problem_IO_3Dfield.cpp
	src/RT_problem_txt_IO.cpp
	src/RT_solver.cpp
	src/Utilities.cpp
	src/RT_utility.cpp
	src/RT_arbitrary_beam.cpp
	src/domain_decomposition_3D.cpp
	src/RT_problem_in.cpp
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
	src/domain_decomposition_3D.hpp
)

### RII LINKAGE
set(RII_INCLUDES "${RII_ROOT_PATH}")
message(STATUS "RII includes: ${RII_INCLUDES}")

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
include_directories(${RII_INCLUDES}/src/riii_approx_two_terms)
include_directories(${RII_INCLUDES}/src/rii_aa_twoterm)
include_directories(${RII_INCLUDES}/src/rii_aa_twoterm_unmag)
if(USE_ACC)
	include_directories(${RII_INCLUDES}/src_mpi/RII_contrib)
	include_directories(${RII_INCLUDES}/src_acc/src)
endif()

### PETSC DEPENDENCY
message(STATUS "Search Petsc")
set(PETSC_TEST_RUNS TRUE)

# Check if PETSC_DIR equals $HOME/petsc
if("$ENV{PETSC_DIR}"
	STREQUAL
	"$ENV{HOME}/petsc"
) # home build in single precision

	set(ENV{PKG_CONFIG_PATH} $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig)

else() # module build in double precision

	set(ENV{PKG_CONFIG_PATH} $ENV{PETSC_DIR}/lib/pkgconfig)

endif()

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
set(RII_SO_BUILD "build" CACHE STRING "Directory name for the librii.so")
set(RII_SO_PATH "${RII_ROOT_PATH}/${RII_SO_BUILD}/librii.so")

set_property(
	TARGET rii
	PROPERTY
		IMPORTED_LOCATION
		"${RII_ROOT_PATH}/${RII_SO_BUILD}/librii.so"
)
message(STATUS "RII path: ${RII_ROOT_PATH}")
message(STATUS "RII SO path ${RII_SO_PATH}")

# import all the dependencies and targets
add_subdirectory(src/GridManager)
target_link_libraries(${EXECUTABLE_NAME} PUBLIC GridManager)
target_link_libraries(${EXECUTABLE_NAME} PUBLIC rii)
target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${Boost_LIBRARIES})
target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${TRIP_MPI_TARGET})
target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${HDF5_LIBRARIES})
target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${HDF5_HL_LIBRARIES})
target_sources(${EXECUTABLE_NAME} PUBLIC ${SOURCES_FILES})
target_include_directories(${EXECUTABLE_NAME}
	PUBLIC
		src
		${CMAKE_CURRENT_BINARY_DIR} # Add binary dir for configure.h
)
target_include_directories(${EXECUTABLE_NAME} SYSTEM
	PUBLIC ${PETSC_INCLUDE_DIRS}
)
target_link_libraries(${EXECUTABLE_NAME} PUBLIC PkgConfig::PETSC)
target_compile_definitions(${EXECUTABLE_NAME} PUBLIC ${PETSC_DEFINITIONS})
# target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${PETSC_LIBRARIES})

message(STATUS "PETSC_INCLUDES = ${PETSC_INCLUDE_DIRS}")

if(NOT yaml-cpp_FOUND)
	target_link_libraries(${EXECUTABLE_NAME} PUBLIC yaml-cpp::yaml-cpp)
else()
	target_link_libraries(${EXECUTABLE_NAME} PUBLIC yaml-cpp)
endif()

# Build thdf/libhdf_output.a using its Makefile
if(TRUE)
	add_subdirectory(src/thdf)
	target_link_libraries(${EXECUTABLE_NAME} PRIVATE hdf_output)
	target_include_directories(${EXECUTABLE_NAME} PUBLIC src/thdf/)

	add_subdirectory(src/atmos_input_TRIP)
	target_link_libraries(${EXECUTABLE_NAME} PRIVATE hdf_atmos)
	target_include_directories(${EXECUTABLE_NAME} PUBLIC src/atmos_input_TRIP/)
endif()

###### make test ######
add_executable(solar_3D_tests ${CMAKE_SOURCE_DIR}/tests/test.cpp)

# import all the dependencies and targets
target_link_libraries(solar_3D_tests
	PUBLIC
		rii
		${Boost_LIBRARIES}
		${TRIP_MPI_TARGET}
)

target_link_libraries(solar_3D_tests PRIVATE hdf_output)
target_include_directories(solar_3D_tests PUBLIC src/thdf/)

target_link_libraries(solar_3D_tests PRIVATE hdf_atmos)
target_include_directories(solar_3D_tests PUBLIC src/atmos_input_TRIP/)

# target_link_libraries(solar_3D_tests PUBLIC ${HDF5_LIBRARIES})
# target_link_libraries(solar_3D_tests PUBLIC ${HDF5_HL_LIBRARIES})

target_sources(solar_3D_tests PUBLIC ${SOURCES_FILES})
target_include_directories(solar_3D_tests
	PUBLIC
		src
		${CMAKE_CURRENT_BINARY_DIR} # Add binary dir for configure.h
)
target_include_directories(solar_3D_tests PUBLIC ${HDF5_INCLUDE_DIRS})

target_include_directories(solar_3D_tests SYSTEM PUBLIC ${PETSC_INCLUDE_DIRS})
target_link_libraries(solar_3D_tests PUBLIC PkgConfig::PETSC)
target_compile_definitions(solar_3D_tests PUBLIC ${PETSC_DEFINITIONS})

if(NOT yaml-cpp_FOUND)
	target_link_libraries(solar_3D_tests PUBLIC yaml-cpp::yaml-cpp)
else()
	target_link_libraries(solar_3D_tests PUBLIC yaml-cpp)
endif()

target_link_libraries(solar_3D_tests PUBLIC GridManager)
target_link_libraries(solar_3D_tests PRIVATE hdf_atmos)

# Include thdf headers
target_include_directories(solar_3D_tests PUBLIC src/thdf/)
target_include_directories(solar_3D_tests PUBLIC src/atmos_input_TRIP)

# Link the static library
target_link_libraries(solar_3D_tests
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/thdf/libhdf_output.a
)

# make does not build test
set_target_properties(solar_3D_tests PROPERTIES EXCLUDE_FROM_ALL TRUE)

# make test build test
add_custom_target(build_test DEPENDS solar_3D_tests)
