# PORTA TRIP Demo

## `write_main_demo_TRIP(int argc, char **argv)`

Serial demo writer that creates an example atmosphere file named `atmosphere.h5`.

Related source files:

- Serial demo implementation: [`atmos_input_TRIP/hdf_atmos_example_trip.c`](./hdf_atmos_example_trip.c)
- MPI demo implementation: [`atmos_input_TRIP/hdf_atmos_mpi.c`](./hdf_atmos_mpi.c)
- Atmosphere API implementation: [`atmos_input_TRIP/hdf_atmos.c`](./hdf_atmos.c)
- Atmosphere API header: [`atmos_input_TRIP/hdf_atmos.h`](./hdf_atmos.h)

### Prototype

```c
int write_main_demo_TRIP(int argc, char **argv);
```

### Parameters

- `argc`: currently unused.
- `argv`: currently unused.

### What It Does

1. Creates or truncates `atmosphere.h5` with `H5Fcreate`.
2. Defines a demo grid with dimensions `N_x = 30`, `N_y = 50`, `N_z = 20`.
3. Builds a height array and writes geometry using `write_geometry_to_hdf5`.
4. Writes a two-level atom description using `write_atom_to_hdf5`.
5. Creates the atmosphere dataset with `THDF_create_atmosphere_dataset`.
6. Populates atmospheric points with `create_atmos_data_point(i, j, k)`.
7. Writes data in `(1, 1, N_z)` blocks via `THDF_write_atmosphere_data_to_dataset`.
8. Closes HDF5 handles and frees allocated memory.

Useful API references:

- HDF5 `H5Fcreate`: <https://portal.hdfgroup.org/documentation/hdf5/latest/group___h5_f.html>
- C `EXIT_SUCCESS` / `EXIT_FAILURE`: <https://en.cppreference.com/w/c/program/EXIT_status>

### How to write 3D parameters in the dataset (Serial version)

Reference implementation: [`atmos_input_TRIP/hdf_atmos_example_trip.c`](./hdf_atmos_example_trip.c).

Serial workflow:

1. Create the dataset handler with `THDF_create_atmosphere_dataset(file, N_x, N_y, N_z)`.
2. Fill a local buffer with row-major indexing `(x, y, z)`, where `z` is the fastest index.
3. Write one block per `(i, j)` with start `(i, j, 0)` and count `(1, 1, N_z)`.
4. Close the dataset handler with `THDF_close_atmosphere_dataset`.

Example write call:

```c
THDF_write_atmosphere_data_to_dataset(
        atmos_dset,
        &atmos_data[idx_base],  // First point of the block (i, j, 0)
        i, j, 0,                // Start indices
        1, 1, N_z);             // Block size
```

Geometry setup example: [`atmos_input_TRIP/hdf_atmos_example_trip.c`](./hdf_atmos_example_trip.c).

See also:

- Atmosphere writer API definitions: [`atmos_input_TRIP/hdf_atmos.h`](./hdf_atmos.h)
- Atmosphere writer API implementation: [`atmos_input_TRIP/hdf_atmos.c`](./hdf_atmos.c)


### Return Value

- `EXIT_SUCCESS` on success.
- `EXIT_FAILURE` if file creation, memory allocation, or HDF5 writes fail.

### Notes

- This function is intended for **serial** execution.
- `main()` checks launcher size and exits if the executable is started with multiple MPI ranks.

### How to write 3D parameters in the dataset (MPI version)

Full MPI example: [`atmos_input_TRIP/hdf_atmos_mpi.c`](./hdf_atmos_mpi.c).

MPI workflow summary:

1. Rank `0` writes shared metadata (atom, frequency grid, geometry).
2. All ranks participate in parallel writes for the 3D atmospheric dataset.
3. If only a subset of ranks should write, create a dedicated communicator with `MPI_Comm_split` and use that communicator consistently for I/O.

Important rule: processes in the writer communicator must enter the HDF5 parallel write calls collectively.

Useful references:

- `MPI_Comm_split` documentation: <https://www.open-mpi.org/doc/current/man3/MPI_Comm_split.3.php>
- Parallel HDF5 overview: <https://portal.hdfgroup.org/documentation/hdf5/latest/_intro_par_h_d_f5.html>
- Parallel HDF5 collective semantics: <https://portal.hdfgroup.org/documentation/hdf5/latest/collective_calls.html>

### Coming soon 
Continuum .....






