# HDF5 Atmosphere Demo Makefile
# Compiler and flags
CC = gcc
CXX = g++
MPICC = mpicc
CFLAGS = -Wall -Wextra -std=c99 -O2 -MMD -MP
CXXFLAGS = -Wall -Wextra -std=c++11 -O2 -MMD -MP
LDFLAGS = -lhdf5 -lhdf5_hl -lm
MPI_LDFLAGS = -lhdf5 -lhdf5_hl -lm

# Directory structure
BUILD_DIR = build

# Project settings
TARGET = hdf_atmos_example
TRIP_TARGET = hdf_atmos_example_trip
CPP_TARGET = hdf_atmos_cpp_example
MPI_TARGET = hdf_atmos_mpi_example
SOURCES = hdf_atmos.c hdf_atmos_functions.c hdf_atmos_continuum.c
TRIP_SOURCES = hdf_atmos_example_trip.c hdf_atmos_functions.c
MPI_SOURCES = hdf_atmos_mpi.c hdf_atmos_functions.c hdf_atmos_continuum.c
CPP_SOURCES = hdf_atmos_cpp_example.cpp
CPP_C_SOURCES = hdf_atmos_functions.c
OBJECTS = $(addprefix $(BUILD_DIR)/,$(SOURCES:.c=.o))
TRIP_OBJECTS = $(addprefix $(BUILD_DIR)/,$(TRIP_SOURCES:.c=.o))
MPI_OBJECTS = $(addprefix $(BUILD_DIR)/,$(MPI_SOURCES:.c=.o))
CPP_OBJECTS = $(addprefix $(BUILD_DIR)/,$(CPP_SOURCES:.cpp=.o))
CPP_C_OBJECTS = $(addprefix $(BUILD_DIR)/,$(CPP_C_SOURCES:.c=.o))
DEPS = $(OBJECTS:.o=.d) $(TRIP_OBJECTS:.o=.d) $(MPI_OBJECTS:.o=.d) $(CPP_OBJECTS:.o=.d)
HEADERS = hdf_atmos.h hdf_atmos_cpp.h

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Build target
$(TARGET): $(BUILD_DIR) $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

# Build trip target
$(TRIP_TARGET): $(BUILD_DIR) $(TRIP_OBJECTS)
	$(CC) $(TRIP_OBJECTS) -o $@ $(LDFLAGS)

# C++ build target
$(CPP_TARGET): $(BUILD_DIR) $(CPP_OBJECTS) $(CPP_C_OBJECTS)
	$(CXX) $(CPP_OBJECTS) $(CPP_C_OBJECTS) -o $@ $(LDFLAGS)

# MPI build target
$(MPI_TARGET): $(BUILD_DIR) $(MPI_OBJECTS) 
	$(MPICC) $(MPI_OBJECTS) -o $@ $(MPI_LDFLAGS)

# Object files
$(BUILD_DIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# C++ object files
$(BUILD_DIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# MPI object files (use MPICC for mpi files)
$(BUILD_DIR)/hdf_atmos_mpi.o: hdf_atmos_mpi.c
	$(MPICC) $(CFLAGS) -c $< -o $@

# Phony targets
.PHONY: all clean run debug help mpi run-mpi cpp run-cpp trip

all: $(TARGET)

# Build C++ version
cpp: $(CPP_TARGET)

# Build trip version
trip: $(TRIP_TARGET)

# Build MPI version
mpi: $(MPI_TARGET)

# Run the program
run: $(TARGET)
	./$(TARGET)

# Run C++ program
run-cpp: $(CPP_TARGET)
	./$(CPP_TARGET)

# Run MPI program (example with 4 processes)
run-mpi: $(MPI_TARGET)
	mpirun -np 4 ./$(MPI_TARGET)

# Debug build
debug: CFLAGS += -g -DDEBUG
debug: CXXFLAGS += -g -DDEBUG
debug: clean $(TARGET)

# Clean build artifacts
clean:
	rm -rfv $(BUILD_DIR) $(TARGET) $(TRIP_TARGET) $(CPP_TARGET) $(MPI_TARGET) atmosphere.h5 atmosphere_mpi.h5 *.a *.h5

# Show help
help:
	@echo "Available targets:"
	@echo "  all     - Build the C project (default)"
	@echo "  trip    - Build Tanausu C version"
	@echo "  cpp     - Build C++ version"
	@echo "  mpi     - Build MPI version"
	@echo "  clean   - Remove build artifacts"
	@echo "  run     - Build and run the C program"
	@echo "  run-cpp - Build and run C++ program"
	@echo "  run-mpi - Build and run MPI program with 4 processes"
	@echo "  debug   - Build with debug symbols"
	@echo "  help    - Show this help message"

# Include dependency files
-include $(DEPS)
