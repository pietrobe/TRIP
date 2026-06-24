COMPILER ?= standard

CC ?= gcc
MPICC ?= mpicc

# Base flags applied to all compilations
CFLAGS += -O3 -Wall -Wextra -march=native -mtune=native -I.

ifeq ($(COMPILER),cray)
	# On Cray systems, 'cc' acts as both the C and MPI compiler,
	# and automatically links loaded modules (like HDF5) seamlessly!
	CC = cc
	MPICC = cc
	CFLAGS += -std=gnu18
	HDF5_CFLAGS =
	HDF5_LIBS =
else ifeq ($(COMPILER),nvidia)
	# NVIDIA HPC SDK specific flags
	CFLAGS += -std=gnu18 -D_Float128=__float128
else
	# Standard compilation flags (GCC, Clang, etc.)
	CFLAGS += -std=c18
endif

# HDF5 configuration
HDF5_CFLAGS = $(shell pkg-config --cflags hdf5 2>/dev/null || h5cc -showconfig 2>/dev/null | grep -i cflags)
HDF5_LIBS = $(shell pkg-config --libs hdf5 hdf5_hl 2>/dev/null || echo "-lhdf5 -lhdf5_hl")

CFLAGS += $(HDF5_CFLAGS)
LDFLAGS = $(HDF5_LIBS) -lm

# Architecture-specific build directory: build/x86_64, build/aarch64, etc.
# This prevents stale cross-architecture objects from causing linker failures
# when the same checkout is built on different machines (e.g. x86 laptop + aarch64 HPC).
ARCH      := $(shell uname -m)
BUILD_DIR  = build/$(ARCH)

# Source files
SRCS = output_utils.c hdf_output_field.c \
	hdf_output_ar_field.c \
	hdf_output_JKQ.c hdf_output_JKQ_CRD.c hdf_output_3D_field.c \
	hdf_output_example_3D.c \
	hdf_fapl_from_env.c \
	hdf_output_domain_decomposition_3D.c \
	hdf_output_zslice_3D.c \
	hdf_output_single_zslice_3D.c \
	hdf_output_example_zslice_3D.c  \
	hdf_output_example_multi_zslice_3D.c \
	hdf_output_example_single_zslice_3D.c

OBJS = $(SRCS:%.c=$(BUILD_DIR)/%.o)
DEPS = $(SRCS:%.c=$(BUILD_DIR)/%.d)

# Library and executables live in the arch-specific build dir too,
# so 'make mpi CC=mpicc' on aarch64 never picks up an x86_64 .a
LIB     = $(BUILD_DIR)/libhdf_output.a
EXE     = hdf_output_example
DEMO    = hdf_output_demo
MPI_EXE = hdf_output_mpi_example

.PHONY: all clean demo lib mpi

all: example

lib: $(LIB)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

example: $(EXE)

$(EXE): hdf_output_example.c $(LIB)
	$(MPICC) $(CFLAGS) -o $@ $< -L$(BUILD_DIR) -lhdf_output $(LDFLAGS)

mpi: CC = $(MPICC)
mpi: $(MPI_EXE)

$(MPI_EXE): hdf_output_mpi_example.c $(LIB)
	$(MPICC) $(CFLAGS) -o $@ $< -L$(BUILD_DIR) -lhdf_output $(LDFLAGS)

demo: $(DEMO)

$(DEMO): hdf_output_demo.c $(LIB) ../atmos/libhdf_atmos.a
	$(MPICC) $(CFLAGS) -o $@ $< -L$(BUILD_DIR) -lhdf_output $(LDFLAGS)

$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rfv $(BUILD_DIR) $(EXE) $(DEMO) $(MPI_EXE) *.h5 *.png 3d_field/*.h5 3d_field

# Include dependency files
-include $(DEPS)
