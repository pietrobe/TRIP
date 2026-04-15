#ifndef __THDF_OUTPUT_FIELD_H__
#define __THDF_OUTPUT_FIELD_H__

#include <float.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <hdf5.h>
#include <hdf5_hl.h>
#include <mpi.h>

#include "../RT_types.hpp"

#ifdef __cplusplus
extern "C" {
#endif

typedef Real THDF_float_t;
typedef Real THDF_n_float_t;

typedef float THDF_float32_t;

typedef enum { HDF_OUT_NONE, HDF_OUT_FLOAT32, HDF_OUT_FLOAT64 } THDF_float_type_t;
typedef enum { HDF_OUT_NORM_NONE, HDF_OUT_NORM_FLOAT32, HDF_OUT_NORM_FLOAT64 } THDF_n_float_type_t;

/**
 * @brief Structure to hold the Stokes parameters field data at a given spatial location and viewing angle
 * that will be written to or read from HDF5 files.
 */
typedef struct {
  int index_i;        // horizontal x index
  int index_j;        // horizontal y index
  int index_incl;     // viewing angle index
  int index_azimuth;  // azimuthal angle index

  // The data must be stored in a ROW MAJOR 5D contiguous array:
  // with memory layout match data layout: [x, y, incl, azimuth, frequencies]
  THDF_float_t *stokes_I;   //
  THDF_float_t *stokes_QI;  //
  THDF_float_t *stokes_UI;  //
  THDF_float_t *stokes_VI;  //

  THDF_n_float_t *norm_multiplier_I;   // TODO: or to evaluate later
  THDF_n_float_t *norm_multiplier_QI;  // Normalization multipliers for each Stokes parameter
  THDF_n_float_t *norm_multiplier_UI;  // At the moment not used
  THDF_n_float_t *norm_multiplier_VI;  //

} THDF_field_t;

/**
 * @brief Structure to hold the frequencies grid information for output fields.
 */
typedef struct {
  int     N_frequencies;
  double *frequencies;  // in Hz a vector of size N_frequencies
} THDF_frequencies_grid_t;

/**
 * @brief Structure to hold the angular grid information for output fields.
 */
typedef struct {
  int     N_directions;
  int     N_inclination_angles;
  int     N_azimuthal_angles;
  double *inclination_angles;  // in radians, size N_directions
  double *azimuthal_angles;    // in radians, size N_directions
  int    *inclinations_indices;
  int    *azimuthal_indices;
} THDF_angular_grid_t;

typedef struct {
  int     N_x;
  int     N_y;
  int     N_z;
  double *heights;
  double  delta;
} THDF_geometry_3D_t;

/**
 * @brief Structure to handle HDF5 output field dataset operations.
 */
typedef struct {
  hid_t file_id;
  hid_t group_id;

  hid_t dataset_id_I;
  hid_t dataset_id_Q;
  hid_t dataset_id_U;
  hid_t dataset_id_V;

  hid_t dataspace_id_I;
  hid_t dataspace_id_Q;
  hid_t dataspace_id_U;
  hid_t dataspace_id_V;

  hid_t datatype_id;
  int   N_x;
  int   N_y;
  int   N_incl;
  int   N_azimuth;
  int   N_frequencies;
  bool  is_open;
} THDF_field_handler_t;

/**
 * @brief Get the floating-point type used for Stokes parameters (e.g., float32 or float64).
 *
 * @return The THDF_float_type_t enumeration value indicating the floating-point type.
 */
THDF_float_type_t
THDF_get_float_type(void);

/**
 * @brief Get the HDF5 datatype identifier corresponding to the floating-point type for Stokes parameters.
 *
 * @return The HDF5 hid_t datatype identifier (e.g., H5T_NATIVE_FLOAT or H5T_NATIVE_DOUBLE).
 */
hid_t
THDF_get_hdf_float_datatype(void);

hid_t
THDF_get_hdf_float32_datatype(void);

/**
 * @brief Get the floating-point type used for normalization multipliers (e.g., float32 or float64).
 *
 * @return The THDF_n_float_type_t enumeration value indicating the normalization floating-point type.
 */
THDF_n_float_type_t
THDF_get_n_float_type(void);

/**
 * @brief Get the HDF5 datatype identifier corresponding to the floating-point type for normalization multipliers.
 *
 * @return The HDF5 hid_t datatype identifier for normalization multipliers.
 */
hid_t
THDF_get_hdf_n_float_datatype(void);

/**
 * @brief Create and initialize an HDF5 field handler for MPI-parallel output.
 *
 * Allocates and initializes a handler structure for writing Stokes parameter fields
 * to an HDF5 file in a parallel (MPI) context. The handler manages datasets and dataspaces
 * for the output fields.
 *
 * @param file           HDF5 file identifier (must be open).
 * @param N_x            Number of grid points in the x direction.
 * @param N_y            Number of grid points in the y direction.
 * @param N_incl         Number of inclination angles.
 * @param N_azimuth      Number of azimuthal angles.
 * @param N_frequencies  Number of frequency points.
 * @return Pointer to the initialized THDF_field_handler_t structure, or NULL on failure.
 */
THDF_field_handler_t *
THDF_create_field_handler_mpi(hid_t file, int N_x, int N_y, int N_incl, int N_azimuth, int N_frequencies);

/**
 * @brief Close and free resources associated with an HDF5 field handler (MPI version).
 *
 * Closes all open HDF5 objects and frees memory associated with the handler.
 *
 * @param output_dset Pointer to the THDF_field_handler_t structure to close and free.
 */
void
THDF_close_field_handler_mpi(THDF_field_handler_t *output_dset);

/**
 * @brief Write a block of Stokes parameter field data to the HDF5 dataset using the handler (MPI version).
 *
 * Writes a specified sub-block of the field data to the corresponding HDF5 datasets,
 * starting at the given indices and covering the specified counts in each dimension.
 *
 * @param output_dset         Pointer to the initialized field handler.
 * @param output_field        Pointer to the field data to write.
 * @param start_i             Starting index in the x direction.
 * @param start_j             Starting index in the y direction.
 * @param start_incl          Starting index for inclination angles.
 * @param start_azimuth       Starting index for azimuthal angles.
 * @param count_i             Number of grid points to write in the x direction.
 * @param count_j             Number of grid points to write in the y direction.
 * @param count_incl          Number of inclination angles to write.
 * @param count_azimuth       Number of azimuthal angles to write.
 * @param count_frequencies   Number of frequency points to write.
 * @return 0 on success, negative value on error.
 */
int
THDF_write_field_dataset_to_hdf5(THDF_field_handler_t *output_dset,    //
                                 THDF_field_t         *output_field,   //
                                 hsize_t               start_i,        //
                                 hsize_t               start_j,        //
                                 hsize_t               start_incl,     //
                                 hsize_t               start_azimuth,  //
                                 hsize_t               count_i,        //
                                 hsize_t               count_j,        //
                                 hsize_t               count_incl,     //
                                 hsize_t               count_azimuth,  //
                                 hsize_t               count_frequencies);           //

/**
 * @brief Write a complete Stokes parameter field to an HDF5 file.
 *
 * Writes the entire contents of a THDF_field_t structure to the specified HDF5 file,
 * creating and populating the relevant datasets for all Stokes parameters and normalization multipliers.
 *
 * @param file           HDF5 file identifier (must be open for writing).
 * @param output_field   Pointer to the field data to write.
 * @param N_x            Number of grid points in the x direction.
 * @param N_y            Number of grid points in the y direction.
 * @param N_incl         Number of inclination angles.
 * @param N_azimuth      Number of azimuthal angles.
 * @param N_frequencies  Number of frequency points.
 * @return 0 on success, negative value on error.
 */
int
THDF_write_field_to_hdf5(hid_t               file,          //
                         const THDF_field_t *output_field,  //
                         int                 N_x,           //
                         int                 N_y,           //
                         int                 N_incl,        //
                         int                 N_azimuth,     //
                         int                 N_frequencies);                //

/**
 * @brief Read a Stokes parameter field from an HDF5 file at a specific grid location and viewing angle.
 *
 * Reads the Stokes parameters and normalization multipliers for the specified spatial and angular indices
 * from the given HDF5 file, and returns a THDF_field_t structure containing the data.
 *
 * @param file           HDF5 file identifier (must be open for reading).
 * @param N_frequencies  Number of frequency points to read.
 * @param x_i            Grid index in the x direction.
 * @param y_i            Grid index in the y direction.
 * @param incl_i         Index for inclination angle.
 * @param azimuth_i      Index for azimuthal angle.
 * @return A THDF_field_t structure containing the read field data.
 */
THDF_field_t
THDF_read_field_from_hdf5(hid_t     file,           //
                          const int N_frequencies,  //
                          const int x_i,            //
                          const int y_i,            //
                          const int incl_i,         //
                          const int azimuth_i);     //

/**
 * @brief Write the angular grid information to an HDF5 file.
 *
 * @param file HDF5 file identifier (must be open for writing)
 * @param angular_grid Pointer to the angular grid structure to write
 * @return int
 */
int                                                                        //
THDF_write_angular_grid_to_hdf5(hid_t                      file,           //
                                const THDF_angular_grid_t *angular_grid);  //

/**
 * @brief Read the angular grid information from an HDF5 file.
 *
 * @param file HDF5 file identifier (must be open for reading)
 * @return THDF_angular_grid_t
 */
THDF_angular_grid_t
THDF_read_angular_grid_from_hdf5(hid_t file);  //

/**
 * @brief Write the frequencies grid information to an HDF5 file.
 *
 * @param file HDF5 file identifier (must be open for writing)
 * @param frequencies_grid Pointer to the frequencies grid structure to write
 * @return int
 */
int
THDF_write_frequencies_grid_to_hdf5(hid_t                          file,               //
                                    const THDF_frequencies_grid_t *frequencies_grid);  //

/**
 * @brief Read the frequencies grid information from an HDF5 file.
 *
 * @param file HDF5 file identifier (must be open for reading)
 * @return THDF_frequencies_grid_t
 */

int
THDF_write_geometry_3D_to_hdf5(hid_t file, const THDF_geometry_3D_t *geometry);  //

/**
 * @brief Read the frequencies grid information from an HDF5 file.
 *
 * @param file HDF5 file identifier (must be open for reading)
 * @return THDF_frequencies_grid_t
 */
THDF_frequencies_grid_t                            //
THDF_read_frequencies_grid_from_hdf5(hid_t file);  //

//============================================================================
// MPI Utility Functions
//============================================================================

/**
 * @brief Open an HDF5 file in MPI parallel mode.
 *
 * Opens an HDF5 file for parallel access using the specified MPI communicator.
 * The file is opened with read/write access.
 *
 * @param filename Path to the HDF5 file to open.
 * @param mpi_comm MPI communicator to use for parallel I/O.
 * @return HDF5 file identifier on success, negative value on failure.
 */
hid_t
THDF_open_file_MPI(const char *filename, MPI_Comm mpi_comm);

/**
 * @brief Create and open an HDF5 file for sequential I/O
 *
 * Creates a new HDF5 file or truncates an existing file with the given filename.
 * This is for non-MPI sequential file creation.
 *
 * @param filename Path to the HDF5 file to create.
 * @return HDF5 file identifier on success, negative value on failure.
 */
hid_t
THDF_open_file(const char *filename);  //

void                             //
THDF_close_file(hid_t file_id);  //

//============================================================================
// Empty structure builders
//============================================================================

THDF_field_t
THDF_make_empty_field(void);  //

THDF_angular_grid_t
THDF_make_empty_angular_grid(void);  //

THDF_frequencies_grid_t
THDF_make_empty_frequencies_grid(void);  //

/* Clear / free helpers: free internal allocation and reset fields */
void
THDF_clear_field(THDF_field_t *output_field);

void
THDF_clear_angular_grid(THDF_angular_grid_t *angular_grid);

void
THDF_clear_frequencies_grid(THDF_frequencies_grid_t *frequencies_grid);

// Macro for HDF5 write error handling
#define HDF5_CHECK_WRITE(status, param_name, memspace)                                \
  do {                                                                                \
    if ((status) < 0) {                                                               \
      fprintf(stderr, "Error writing %s data to output field dataset\n", param_name); \
      H5Sclose(memspace);                                                             \
      return -1;                                                                      \
    }                                                                                 \
  } while (0)

// Macro for HDF5 write error handling
#define HDF5_CHECK_WRITE_OF(status, param_name, memspace, return_value)               \
  do {                                                                                \
    if ((status) < 0) {                                                               \
      fprintf(stderr, "Error writing %s data to output field dataset\n", param_name); \
      H5Sclose(memspace);                                                             \
      return return_value;                                                            \
    }                                                                                 \
  } while (0)

#ifdef __cplusplus
}
#endif

#endif  // __THDF_OUTPUT_FIELD_H__