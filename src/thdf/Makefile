# Compiler and flags
CC ?= gcc
MPICC ?= mpicc
CFLAGS += -Wall -Wextra -O3 -std=c99 --pedantic -I.

# HDF5 configuration (use pkg-config or hdf5 compiler wrapper)
HDF5_CFLAGS = $(shell pkg-config --cflags hdf5 2>/dev/null || h5cc -showconfig | grep -i cflags)
HDF5_LIBS = $(shell pkg-config --libs hdf5 hdf5_hl 2>/dev/null || echo "-lhdf5 -lhdf5_hl")

CFLAGS += $(HDF5_CFLAGS)
LDFLAGS = $(HDF5_LIBS) -lm

# Build directory
BUILD_DIR = build

# Source files
SRCS = hdf_output_field.c hdf_output_ar_field.c
OBJS = $(SRCS:%.c=$(BUILD_DIR)/%.o)
DEPS = $(SRCS:%.c=$(BUILD_DIR)/%.d)

# Library target
LIB = libhdf_output.a
EXE = hdf_output_example

# Demo executable
DEMO = hdf_output_demo

# MPI executable
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
	$(CC) $(CFLAGS) -o $@ $< -L. -lhdf_output $(LDFLAGS)

mpi: $(MPI_EXE)
$(MPI_EXE): hdf_output_mpi_example.c $(LIB)
	$(MPICC) $(CFLAGS) -o $@ $< -L. -lhdf_output $(LDFLAGS)
demo: $(DEMO)

$(DEMO): hdf_output_demo.c $(LIB) ../atmos/libhdf_atmos.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lhdf_output $(LDFLAGS)
$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rfv  $(LIB) $(EXE) $(DEMO) $(MPI_EXE) *.h5 $(BUILD_DIR)

# Include dependency files
-include $(DEPS)