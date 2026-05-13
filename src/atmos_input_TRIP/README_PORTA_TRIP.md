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
6. Populates atmospheric points with `create_atmos_data_point(i, j, k, &atom)`.
7. Writes data in `(1, 1, N_z)` blocks via `THDF_write_atmosphere_data_to_dataset`.
8. Closes HDF5 handles and frees allocated memory.

Useful API references:

- HDF5 `H5Fcreate`: <https://portal.hdfgroup.org/documentation/hdf5/latest/group___h5_f.html>
- C `EXIT_SUCCESS` / `EXIT_FAILURE`: <https://en.cppreference.com/w/c/program/EXIT_status>

### Return Value

- `EXIT_SUCCESS` on success.
- `EXIT_FAILURE` if file creation, memory allocation, or HDF5 writes fail.

### Notes

- This function is intended for **serial** execution.
- `main()` checks launcher size and exits if the executable is started with multiple MPI ranks.

### How to Write 3D Parameters (Serial)

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

See also:

- Atmosphere writer API definitions: [`atmos_input_TRIP/hdf_atmos.h`](./hdf_atmos.h)
- Atmosphere writer API implementation: [`atmos_input_TRIP/hdf_atmos.c`](./hdf_atmos.c)

## How to Write 3D Parameters (MPI)

Full MPI example: [`atmos_input_TRIP/hdf_atmos_mpi.c`](./hdf_atmos_mpi.c).

The parallel workflow has three distinct phases:

1. **All ranks** create the file collectively with a parallel property list, then immediately close it.
2. **Rank 0 only** reopens the file with serial I/O and writes all shared metadata (geometry, atom).
3. **All ranks** reopen the file with parallel I/O, each rank writes its own spatial slice.

### Phase 1 — Create the file collectively

Every rank must participate in `H5Fcreate`. The file-access property list must be
configured for MPI-IO before the call:

```c
hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

hid_t file = H5Fcreate("atmosphere_mpi.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
H5Pclose(plist_id);

// Close immediately — rank 0 will reopen with serial I/O for metadata
H5Fclose(file);
MPI_Barrier(MPI_COMM_WORLD);
```

### Phase 2 — Rank 0 writes metadata (serial)

Only rank 0 writes shared metadata. It opens the file with the default (serial)
property list and calls the metadata writers:

```c
if (mpi_rank == 0) {
    file = H5Fopen("atmosphere_mpi.h5", H5F_ACC_RDWR, H5P_DEFAULT);

    write_geometry_to_hdf5(file, N_x, N_y, N_z, heights);
    write_atom_to_hdf5(file, &atom);

    H5Fclose(file);
}
MPI_Barrier(MPI_COMM_WORLD);
```

All other ranks wait at the barrier before proceeding to the parallel write.

### Phase 3 — All ranks write in parallel

#### Domain decomposition

The example divides the `x` axis evenly across ranks, distributing any remainder
to the first ranks:

```c
int local_N_x = N_x / mpi_size;
int remainder = N_x % mpi_size;
int start_i   = mpi_rank * local_N_x + (mpi_rank < remainder ? mpi_rank : remainder);
if (mpi_rank < remainder)
    local_N_x++;
```

Each rank owns `local_N_x` slices starting at global index `start_i`.

#### Reopen the file for parallel I/O

```c
plist_id = H5Pcreate(H5P_FILE_ACCESS);
H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
file = H5Fopen("atmosphere_mpi.h5", H5F_ACC_RDWR, plist_id);
H5Pclose(plist_id);
```

#### Create the MPI dataset handle (collective)

All ranks call `create_atmosphere_dataset_mpi` collectively:

```c
THDF_atmos_dataset_handler_t *atmos_dset =
    create_atmosphere_dataset_mpi(file, N_x, N_y, N_z, H5P_DEFAULT);
```

This is a collective operation — every rank in the communicator must call it.

#### Fill local atmosphere data

Each rank fills its own `local_atmos_data` buffer with global indices:

```c
for (int i = 0; i < local_N_x; i++) {
    for (int j = 0; j < N_y; j++) {
        for (int k = 0; k < N_z; k++) {
            int local_idx   = i * N_y * N_z + j * N_z + k;
            int global_i    = start_i + i;
            local_atmos_data[local_idx] =
                create_atmos_data_point_mpi(global_i, j, k, mpi_rank, &atom);
        }
    }
}
```

The `index_i`, `index_j`, `index_k` fields inside each `THDF_atmos_t` must
hold the **global** grid coordinates, not local offsets.

#### Set the data-transfer property list

```c
hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
```

Use `H5FD_MPIO_COLLECTIVE` if all ranks always write (often faster due to MPI-IO
optimisations); use `H5FD_MPIO_INDEPENDENT` when some ranks may skip a write.

#### Write the local atmosphere slice

```c
write_atmosphere_data_to_dataset_mpi(
    atmos_dset,
    local_atmos_data,
    start_i, 0, 0,          // global start (i, j, k)
    local_N_x, N_y, N_z,   // count
    plist_id);
```

#### Clean up

```c
close_atmosphere_dataset_mpi(atmos_dset);
H5Pclose(plist_id);
H5Fclose(file);
```

### Verification (rank 0 only)

After `MPI_Finalize` is not yet called, rank 0 can read back the file to
validate it:

```c
if (mpi_rank == 0) {
    THDF_read_main_demo_from_file("atmosphere_mpi.h5", false);
}
```

### Important rules

- Every rank in the communicator must call `create_atmosphere_dataset_mpi` and `write_atmosphere_data_to_dataset_mpi` — even ranks with zero local points (use `count_i = 0` in that case).
- Do **not** mix serial and parallel file handles for the same file simultaneously across ranks.
- The `THDF_atmos_t` grid indices (`index_i`, `index_j`, `index_k`) must be the **global** coordinates, not local buffer offsets.
- If only a subset of ranks needs to write, create a dedicated communicator with `MPI_Comm_split` and pass it to `H5Pset_fapl_mpio`; all other ranks must not call the collective HDF5 functions on that handle.

Useful references:

- `MPI_Comm_split` documentation: <https://www.open-mpi.org/doc/current/man3/MPI_Comm_split.3.php>
- Parallel HDF5 overview: <https://portal.hdfgroup.org/documentation/hdf5/latest/_intro_par_h_d_f5.html>
- Parallel HDF5 collective semantics: <https://portal.hdfgroup.org/documentation/hdf5/latest/collective_calls.html>

## Continuum

Continuum: The current version of THDF supports frequency-dependent continuum opacity, but this feature is not yet used.

## `THDF_atmos_t`

Main per-grid-point structure for atmosphere data, defined in [`atmos_input_TRIP/hdf_atmos.h`](./hdf_atmos.h).

Each element of the 3D dataset corresponds to one `THDF_atmos_t` value stored at grid position `(index_i, index_j, index_k)`.

### Definition

```c
typedef struct {
  int16_t index_i;
  int16_t index_j;
  int16_t index_k;

  float  temperature;
  double nH;
  double vmicro;
  double Doppler_width;
  double damping;
  double pop_lower_level;
  double pop_upper_level;
  double Cul;
  double Qel;
  double D2;

  float magnetic_field_x;
  float magnetic_field_y;
  float magnetic_field_z;

  float bulk_velocity_x;
  float bulk_velocity_y;
  float bulk_velocity_z;

  double c_scat_opacity_sigma_c;
  double c_therm_opacity_k_c;
  double c_therm_emissivity_epsilon_cth;
} THDF_atmos_t;
```

### Fields

| Field | Type | Unit | Description |
| --- | --- | --- | --- |
| `index_i` | `int16_t` | — | Grid index along x; must be unique per point |
| `index_j` | `int16_t` | — | Grid index along y; must be unique per point |
| `index_k` | `int16_t` | — | Grid index along z (height); must be unique per point |
| `temperature` | `float` | K | Local temperature |
| `nH` | `double` | cm⁻³ | Hydrogen atom number density |
| `vmicro` | `double` | cm/s | Microturbulent velocity |
| `Doppler_width` | `double` | Hz | Doppler line width |
| `damping` | `double` | — | Damping parameter (dimensionless) |
| `pop_lower_level` | `double` | cm⁻³ | Population of the lower energy level |
| `pop_upper_level` | `double` | cm⁻³ | Population of the upper energy level |
| `Cul` | `double` | s⁻¹ | Rate of inelastic collisions |
| `Qel` | `double` | s⁻¹ | Rate of elastic collisions with neutral perturbers |
| `D2` | `double` | — | D² depolarization rate. Derouich, M. (2019), *The Astrophysical Journal*, 880, 10 (6 pp.), 2019 July 20 |
| `magnetic_field_x` | `float` | G | Magnetic field x component |
| `magnetic_field_y` | `float` | G | Magnetic field y component |
| `magnetic_field_z` | `float` | G | Magnetic field z component |
| `bulk_velocity_x` | `float` | cm/s | Bulk velocity x component |
| `bulk_velocity_y` | `float` | cm/s | Bulk velocity y component |
| `bulk_velocity_z` | `float` | cm/s | Bulk velocity z component |
| `c_scat_opacity_sigma_c` | `double` | — | Continuum scattering opacity σ_c |
| `c_therm_opacity_k_c` | `double` | — | Continuum thermal opacity k_c |
| `c_therm_emissivity_epsilon_cth` | `double` | — | Continuum thermal emissivity ε_cth |

### Continuum fields

The three continuum fields hold values at the **central line frequency**. The total continuum opacity is:

```text
K_c = c_therm_opacity_k_c + c_scat_opacity_sigma_c
```

For blends or frequency-resolved calculations, per-frequency arrays would be needed instead (see [Continuum](#continuum)).
