#ifndef HDF_ATMOS_H
#define HDF_ATMOS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Structure representing a single point in the 3D atmosphere grid
 */
typedef struct {
  int16_t index_i;  // horizontal x index
  int16_t index_j;  // horizontal y index
  int16_t index_k;  // height index   // These must be unique for each grid point

  float  temperature;      // Temperature in K
  double nH;               // Hydrogen atoms density
  double vmicro;           // Microturbulent velocity in cm/s (set to 0)
  double Doppler_width;    // Doppler W Hz
  double damping;          // Damping parameter (dimensionless)
  double pop_lower_level;  // Population of lower energy level (in cm^-3)
  double pop_upper_level;  // Population of upper energy level (in cm^-3)
  double Cul;              // Rate of inelastic collisions (in s^-1)
  double Qel;              // Rate of elastic collisions with neutral perturber (in s^-1)
  double D2;  // D^2 Depolarization Rate: From: Derouich (2019) - The Astrophysical Journal, 880:10 (6pp), 2019 July 20

  float magnetic_field_x;  // in Gauss
  float magnetic_field_y;  // in Gauss
  float magnetic_field_z;  // in Gauss

  float bulk_velocity_x;  // in cm/s
  float bulk_velocity_y;  // in cm/s
  float bulk_velocity_z;  // in cm/s

  // Continuum at central line frequency
  // In case of blends or more accurate calculations,
  // these quantities must be provided at each frequency point
  // with array of size N_frequencies

  // ATTENTION !!!!!!!!!!!!!!!!!!!!!!!
  // the user must calculate /////////
  // Kc = c_therm_opacity_k_c + c_scat_opacity_sigma_c
  double c_scat_opacity_sigma_c;        // sigma_c
  double c_therm_opacity_k_c;           // k_c (new) / Corresponds to k_c (old) - sigma_c
  double c_therm_emissivity_epsilon_c;  // epsilon_cth

} THDF_atmos_t;

/**
 * @brief Structure representing a two-level atom for HDF5 storage
 */
typedef struct {
  int    atomic_number;  // Atomic number
  double atomic_mass;    // Relative atomic mass in amu
  double E_lower;        // in cm^-1
  double E_upper;        // in cm^-1
  double g_lower;        //
  double g_upper;        //
  int    jl2;            // Multiplied by 2
  int    ju2;            //
  double Aul;            //

  double a_coef_D2;  // D^2 = a * nH * (T/5000)^b // Depolarization Rate.
  double b_coef_D2;  // From: Derouich (2019) - The Astrophysical Journal, 880:10 (6pp), 2019 July 20
} THDF_atom_two_levels_t;

/**
 * @brief Maximum number of atomic levels for placeholder arrays in THDF_atom_two_terms_t
 */
#define MAX_ATOM_TERMS_PLACEHOLDER 16

/**
 * @brief Structure representing a two-term atom for HDF5 storage
 */
typedef struct {
  int    atomic_number;  // Atomic number
  double atomic_mass;    // Relative atomic mass in amu
  double Aul;            // Einstein A coefficient for spontaneous emission (in s^-1)
  int    L2_upper;       // Upper level orbital angular momentum (multiplied by 2)
  int    L2_lower;       // Lower level orbital angular momentum (multiplied by 2)
  int    S2;             // total spin multiplied by 2

  int    E_size;                                // Number of energy levels (should not exceed MAX_ATOM_TERMS_PLACEHOLDER)
  double E_lower[MAX_ATOM_TERMS_PLACEHOLDER];  // Lower level energy in cm^-1
  double E_upper[MAX_ATOM_TERMS_PLACEHOLDER];  // Upper level energy in cm^-1
} THDF_atom_two_terms_t;

/**
 * @brief Structure representing the frequency grid for HDF5 storage
 */
typedef struct {
  int     N_frequencies;
  double *frequency;  // in Hz a vector of size N_frequencies
  int     nu0_index;  // index of central line frequency in frequency array (or the closest one)
  // Remove nu0
} THDF_frequency_grid_t;

typedef struct {
  int     N_z;
  int     N_x;
  int     N_y;
  double *height;  // in cm a vector of size N_height
  double  delta;   // in cm
} THDF_atmos_geometry_t;

typedef struct {
  // EXPERIMENTAL: This struct is a work in progress and at the moment is not used in the example code.
  int index_i;  // horizontal x index
  int index_j;  // horizontal y index
  int index_k;  // height index   // These must be unique for each grid point

  // Normalization parameters
  // 3D mtrices (i,j,k)
  double min_c_scat_opacity_sigma_c;     // min sigma_c
  double max_Vm_c_scat_opacity_sigma_c;  // max Vm sigma_c

  double min_c_therm_opacity_k_c;     // min epsilon_c
  double max_Vm_c_therm_opacity_k_c;  // max Vm epsilon_c

  double min_c_therm_emissivity_epsilon_c;     // min eta_thermal (old: K_c)
  double max_Vm_c_therm_emissivity_epsilon_c;  // max Vm eta_thermal (old: K_c)

  // // arrays storing the continuum data (N_freq)
  // uint16_t *c_scat_opacity_sigma_c;        // sigma_c
  // uint16_t *c_therm_opacity_k_c;  // epsilon_c
  // uint16_t *c_therm_emissivity_epsilon_c;                 // eta_thermal (old: K_c)
} THDF_continuum_t;

typedef struct {
  hid_t file_id;
  hid_t group_id;
  hid_t dataset_id;
  hid_t dataspace_id;
  hid_t datatype_id;
  int   N_x;
  int   N_y;
  int   N_z;
  int   is_open;
} THDF_atmos_dataset_handler_t;

/**
 * @brief Write geometry data to HDF5 file
 *
 * Creates a /geometry group and writes grid dimensions and height scale.
 *
 * @param file HDF5 file handle
 * @param N_x Number of grid points in x direction
 * @param N_y Number of grid points in y direction
 * @param N_z Number of height levels
 * @param heights Array of height values in cm (size N_z)
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int
write_geometry_to_hdf5(hid_t   file,      //
                       int     N_x,       //
                       int     N_y,       //
                       int     N_z,       //
                       double *heights);  //

/**
 * @brief Read geometry data from HDF5 file
 *
 * Reads grid dimensions and height scale from /geometry group.
 * Allocates memory for heights array (caller must free).
 *
 * @param file HDF5 file handle
 * @param N_x Pointer to store number of grid points in x direction
 * @param N_y Pointer to store number of grid points in y direction
 * @param N_z Pointer to store number of height levels
 * @param heights Pointer to store allocated height array (in cm)
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int
THDF_read_geometry_from_hdf5(hid_t    file,  //
                             int     *N_x,   //
                             int     *N_y,   //
                             int     *N_z,   //
                             double **heights);

/**
 * @brief Write atomic transition parameters to HDF5 file
 *
 * Creates /atom group and writes atomic data as individual datasets.
 *
 * @param file HDF5 file handle
 * @param atom Pointer to atom structure with transition parameters
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int
write_atom_to_hdf5(hid_t                   file,   //
                   THDF_atom_two_levels_t *atom);  //

/**
 * @brief Write two-term atom data to HDF5 file
 *
 * Creates /atom_two_terms_data dataset and writes the atom structure to it.
 * atom->E_size must not exceed MAX_ATOM_TERMS_PLACEHOLDER. Entries of
 * E_lower/E_upper beyond E_size are zeroed out before writing.
 *
 * @param file HDF5 file handle
 * @param atom Pointer to atom structure with two-term transition parameters
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int
write_atom_two_terms_to_hdf5(hid_t                  file,   //
                             THDF_atom_two_terms_t *atom);  //

/**
 * @brief Check if the atom data in the HDF5 file is a two-term atom
 */
bool
THDF_is_two_term_atom(hid_t file);  //

/**
 * @brief Read atomic transition parameters from HDF5 file
 *
 * Reads atomic data from /atom group.
 *
 * @param file HDF5 file handle
 * @param atom Pointer to atom structure to populate
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int
THDF_read_atom_from_hdf5(hid_t                   file,  //
                         THDF_atom_two_levels_t *atom);

/**
 * @brief Read two-term atom data from HDF5 file
 *
 * Reads atomic data from /atom_two_terms_data dataset.
 *
 * @param file HDF5 file handle
 * @param atom Pointer to atom structure to populate
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int
THDF_read_atom_two_terms_from_hdf5(hid_t                  file,  //
                                   THDF_atom_two_terms_t *atom);

/**
 * @brief Write frequency grid to HDF5 file
 *
 * Creates /frequency_grid group and writes frequency array.
 *
 * @param file HDF5 file handle
 * @param N_frequencies Number of frequency points
 * @param frequencies Array of frequency values in Hz (size N_frequencies)
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int
write_frequency_grid_to_hdf5(hid_t   file,           //
                             int     N_frequencies,  //
                             double *frequencies);

/**
 * @brief Read frequency grid from HDF5 file
 *
 * Reads frequency array from /frequency_grid group.
 * Allocates memory for frequencies array (caller must free).
 *
 * @param file HDF5 file handle
 * @param N_frequencies Pointer to store number of frequency points
 * @param frequencies Pointer to store allocated frequency array (in Hz)
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int
THDF_read_frequency_grid_from_hdf5(hid_t    file,           //
                                   int     *N_frequencies,  //
                                   double **frequencies);

int
THDF_read_atmos_3D_point_hdf5(hid_t file, THDF_atmos_t *atmos_data,    //
                              int start_i, int start_j, int start_k);  //

/**
 * @brief Write 3D atmosphere data to HDF5 file
 *
 * Writes entire 3D atmosphere dataset at root level.
 * Uses HDF5 compound datatype for hdf_atmos structure.
 *
 * @param file HDF5 file handle
 * @param atmos_data Array of atmosphere data points (size N_x * N_y * N_z)
 * @param N_x Number of grid points in x direction
 * @param N_y Number of grid points in y direction
 * @param N_z Number of height levels
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int
write_atmos_3D_data_to_hdf5(hid_t         file,        //
                            THDF_atmos_t *atmos_data,  //
                            int           N_x,         //
                            int           N_y,         //
                            int           N_z);

/**
 * @brief Read 3D atmosphere data from HDF5 file
 *
 * Reads entire 3D atmosphere dataset from root level.
 * Allocates memory for atmos_data array (caller must free).
 * Reads dimensions from dataspace.
 *
 * @param file HDF5 file handle
 * @param atmos_data Pointer to store allocated atmosphere data array
 * @param N_x Pointer to store number of grid points in x direction
 * @param N_y Pointer to store number of grid points in y direction
 * @param N_z Pointer to store number of height levels
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int
THDF_read_atmos_3D_data_from_hdf5(hid_t          file,        //
                                  THDF_atmos_t **atmos_data,  //
                                  int           *N_x,         //
                                  int           *N_y,         //
                                  int           *N_z);

/**
 * @brief Create HDF5 compound datatype for atom two-level structure
 *
 * Defines HDF5 compound type matching the THDF_atom_two_levels_t C struct.
 * Caller must close the returned type with H5Tclose().
 *
 * @return HDF5 datatype handle, negative on error
 */
hid_t
create_atom_compound_type(void);

/**
 * @brief Create HDF5 compound datatype for two-term atom structure
 *
 * Defines HDF5 compound type matching the THDF_atom_two_terms_t C struct.
 * E_lower and E_upper are stored as fixed-size arrays of length
 * MAX_ATOM_TERMS_PLACEHOLDER. Caller must close the returned type with H5Tclose().
 *
 * @return HDF5 datatype handle, negative on error
 */
hid_t
create_atom_two_terms_compound_type(void);

/**
 * @brief Create HDF5 compound datatype for hdf_atmos structure
 *
 * Creates an HDF5 compound datatype that maps to the hdf_atmos C structure.
 * This datatype is used for reading and writing atmosphere data to HDF5 files.
 *
 * The compound type includes the following fields:
 * - index_i, index_j, index_k: Grid position indices (int16_t)
 * - temperature: Temperature in K (float)
 * - vmicro: Microturbulent velocity in cm/s (double)
 * - damping: Damping parameter, dimensionless (double)
 * - pop_lower_level, pop_upper_level: Level populations in cm^-3 (double)
 * - Cul: Inelastic collision rate in s^-1 (double)
 * - Qel: Elastic collision rate in s^-1 (double)

 * - magnetic_field_x, magnetic_field_y, magnetic_field_z: In Gauss (float)
 * - bulk_velocity_x, bulk_velocity_y, bulk_velocity_z: In cm/s (float)
 * - c_scat_opacity_sigma_c: Continuum scattering opacity (double)
 * - c_therm_opacity_k_c: Continuum thermal emissivity (double)
 * - c_therm_emissivity_epsilon_c: Continuum thermal opacity (double) (old: K_c)
 *
 * @return HDF5 datatype handle (hid_t) on success, negative value on error.
 *         Caller is responsible for closing the type with H5Tclose().
 *
 * @note This function is used internally by read/write functions but can also
 *       be called directly for custom HDF5 operations.
 *
 * @see hdf_atmos
 * @see read_atmos_3D_data_from_hdf5
 * @see write_atmos_3D_data_to_hdf5
 */
hid_t
create_hdf_atmos_compound_type(void);

/**
 * @brief Create atmosphere dataset for partial/incremental writing
 *
 * Creates atmosphere dataset at root with specified dimensions.
 * Returns handle for subsequent partial write operations.
 * Must be closed with close_atmosphere_dataset().
 *
 * @param file HDF5 file handle
 * @param N_x Number of grid points in x direction
 * @param N_y Number of grid points in y direction
 * @param N_z Number of height levels
 * @return Pointer to dataset handle, NULL on error
 */
THDF_atmos_dataset_handler_t *
THDF_create_atmosphere_dataset(hid_t file,  //
                               int   N_x,   //
                               int   N_y,   //
                               int   N_z);

/**
 * @brief Write partial atmosphere data to dataset using hyperslab
 *
 * Writes a rectangular subset of atmosphere data to the dataset.
 * Allows incremental writing of large datasets.
 *
 * @param atmos_dset Dataset handle from create_atmosphere_dataset()
 * @param atmos_data Array of atmosphere data to write
// ...existing code...
 * @param start_j Starting index in y direction
 * @param start_k Starting index in z direction (height)
 * @param count_i Number of points to write in x direction
 * @param count_j Number of points to write in y direction
 * @param count_k Number of points to write in z direction
 * @return 0 on success, -1 on error
 */
int
THDF_write_atmosphere_data_to_dataset(THDF_atmos_dataset_handler_t *atmos_dset,  //
                                      THDF_atmos_t                 *atmos_data,  //
                                      hsize_t                       start_i,     //
                                      hsize_t                       start_j,     //
                                      hsize_t                       start_k,     //
                                      hsize_t                       count_i,     //
                                      hsize_t                       count_j,     //
                                      hsize_t                       count_k);

/**
 * @brief Close atmosphere dataset and release resources
 *
 * Closes all HDF5 handles and frees dataset structure.
 *
 * @param atmos_dset Dataset handle to close
 * @return 0 on success, -1 on error
 */
int
THDF_close_atmosphere_dataset(THDF_atmos_dataset_handler_t *atmos_dset);

/**
 * @brief Create HDF5 compound datatype for hdf_atmos structure
 *
 * Defines HDF5 compound type matching the hdf_atmos C struct.
 * Caller must close the returned type with H5Tclose().
 *
 * @return HDF5 datatype handle, negative on error
 */
hid_t
create_hdf_atmos_compound_type(void);

/**
 * @brief Create and populate atmosphere data point with sample values
 *
 * Generates atmosphere data with physics-based sample calculations
 * for testing and demonstration purposes.
 *
 * @param i Grid index in x direction
 * @param j Grid index in y direction
 * @param k Height level index
 * @return Populated hdf_atmos structure
 */
THDF_atmos_t
create_atmos_data_point(int                           i,  //
                        int                           j,  //
                        int                           k,  //
                        const THDF_atom_two_levels_t *atom);

/**
 * @brief Print all fields of an atmosphere data point
 *
 * Prints a human-readable summary of a THDF_atmos_t structure.
 *
 * @param point Pointer to atmosphere data point to print
 */
void
THDF_print_atmos_point(const THDF_atmos_t *point);

/**
 * @brief Create HDF5 compound datatype for atmosphere data
 *
 * Defines HDF5 compound type matching the hdf_atmos C struct.
 *
 * @return hid_t
 */
hid_t
create_atmos_compound_type();

//============================================================================
// MPI-enabled functions for parallel I/O
//============================================================================

/**
 * @brief Create atmosphere dataset for MPI parallel writing
 *
 * MPI-enabled version that creates dataset with parallel property list.
 * All MPI processes must call this collectively.
 *
 * @param file HDF5 file handle (opened with MPI communicator)
 * @param N_x Number of grid points in x direction
 * @param N_y Number of grid points in y direction
 * @param N_z Number of height levels
 * @param plist_id Parallel property list for collective operations
 * @return Pointer to dataset handle, NULL on error
 */
THDF_atmos_dataset_handler_t *
create_atmosphere_dataset_mpi(hid_t file,  //
                              int   N_x,   //
                              int   N_y,   //
                              int   N_z,   //
                              hid_t plist_id);
// ...existing code...
int
write_atmosphere_data_to_dataset_mpi(THDF_atmos_dataset_handler_t *atmos_dset,  //
                                     THDF_atmos_t                 *atmos_data,  //
                                     hsize_t                       start_i,     //
                                     hsize_t                       start_j,     //
                                     hsize_t                       start_k,     //
                                     hsize_t                       count_i,     //
                                     hsize_t                       count_j,     //
                                     hsize_t                       count_k,     //
                                     hid_t                         plist_id);

/**
 * @brief Close MPI atmosphere dataset
 *
 * MPI-enabled version for closing parallel dataset.
 *
 * @param atmos_dset Dataset handle to close
 * @return 0 on success, -1 on error
 */
int
close_atmosphere_dataset_mpi(THDF_atmos_dataset_handler_t *atmos_dset);

/**
 * @brief Create HDF5 compound datatype (MPI version)
 *
 * MPI-enabled version of compound type creation.
 *
 * @return HDF5 datatype handle, negative on error
 */
hid_t
create_hdf_atmos_compound_type_mpi(void);

/**
 * @brief Create atmosphere data point with MPI rank information
 *
 * MPI-enabled version that includes rank information in generated data.
 * Useful for verifying parallel decomposition.
 *
 * @param i Grid index in x direction
 * @param j Grid index in y direction
 * @param k Height level index
 * @param mpi_rank MPI process rank
 * @return Populated hdf_atmos structure
 */
THDF_atmos_t
create_atmos_data_point_mpi(int                           i,         //
                            int                           j,         //
                            int                           k,         //
                            int                           mpi_rank,  //
                            const THDF_atom_two_levels_t *atom);

//============================================================================
// Demo/test functions
//============================================================================

/**
 * @brief Read and validate HDF5 atmosphere file
 *
 * Comprehensive test function that reads all data from file,
 * validates indices and structure, and optionally prints details.
 *
 * @param filename Path to HDF5 file to read
 * @param verbose If true, print detailed information about file contents
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int
THDF_read_main_demo_from_file(const char *filename,  //
                              const bool  verbose);

#ifdef __cplusplus
}
#endif

#endif  // HDF_ATMOS_H