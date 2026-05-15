#ifndef HDF_ATMOS_CONTINUUM_H
#define HDF_ATMOS_CONTINUUM_H

#include <stdint.h>

#include "hdf_atmos.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Create HDF5 compound datatype for THDF_continuum_t
 *
 * Creates a compound type covering the scalar normalization fields.
 * Spectral arrays are stored in dedicated 4D datasets.
 * Caller must close the returned type with H5Tclose().
 *
 * @param N_frequencies Number of frequency points (reserved for future use)
 * @return HDF5 datatype handle on success, H5I_INVALID_HID on error
 */
hid_t
create_hdf_continuum_compound_type(int N_frequencies);

/**
 * @brief Create and initialise the /continuum group and its datasets
 *
 * Creates (or opens) /continuum group, commits the compound type,
 * writes dimension attributes, and creates the 4D uint16 datasets
 * and the 3D compound dataset.
 *
 * @param file         HDF5 file handle
 * @param N_x          Number of x grid points
 * @param N_y          Number of y grid points
 * @param N_z          Number of z grid points
 * @param N_frequencies Number of frequency points
 * @return 0 on success, -1 on error
 */
int
make_continuum_hdf5_mpi(hid_t file, int N_x, int N_y, int N_z, int N_frequencies);

/**
 * @brief Create 4D uint16 continuum datasets under /continuum
 *
 * Creates (or opens and validates) the following datasets with
 * dimensions [N_x, N_y, N_z, N_frequencies]:
 *   - c_scat_opacity_sigma_c
 *   - c_therm_opacity_k_c
 *   - c_therm_emissivity_epsilon_c (old: c_tot_opacity_K_c)
 *
 * Also writes N_x, N_y, N_z, N_frequencies as attributes.
 *
 * @param file         HDF5 file handle
 * @param N_x          Number of x points
 * @param N_y          Number of y points
 * @param N_z          Number of z points
 * @param N_frequencies Number of frequency points
 * @return 0 on success, -1 on error
 */
int
THDF_create_continuum_4D_matrices(hid_t file, int N_x, int N_y, int N_z, int N_frequencies);

/**
 * @brief Normalize continuum spectra to uint16 and compute min/range scalars
 *
 * Computes per-quantity min and (max - min), normalises each frequency
 * sample to [0, UINT16_MAX], and returns a THDF_continuum_t carrying
 * the normalization parameters.
 *
 * @param sigma_c      Raw scattering opacity array       (size N_frequencies)
 * @param epsilon_c    Raw thermal emissivity array        (size N_frequencies)
 * @param c_therm_emissivity_epsilon_c Raw thermal opacity array          (size N_frequencies)
 * @param N_frequencies Number of frequency points
 * @param norm_sigma_c       Output uint16 scattering opacity      (size N_frequencies)
 * @param norm_epsilon_c     Output uint16 thermal emissivity      (size N_frequencies)
 * @param norm_c_therm_emissivity_epsilon_c Output uint16 thermal opacity         (size N_frequencies)
 * @return THDF_continuum_t carrying min/range normalization scalars
 */
THDF_continuum_t
normalize_continuum_data(const double *sigma_c, const double *epsilon_c, const double *c_therm_emissivity_epsilon_c,
                         int N_frequencies, uint16_t *norm_sigma_c, uint16_t *norm_epsilon_c,
                         uint16_t *norm_c_therm_emissivity_epsilon_c);

/**
 * @brief Reconstruct double-precision spectra from uint16 and normalization scalars
 *
 * Inverts the normalization applied by normalize_continuum_data().
 *
 * @param cont           Normalization scalars from normalize_continuum_data()
 * @param norm_sigma_c       Normalized scattering opacity       (size N_frequencies)
 * @param norm_epsilon_c     Normalized thermal emissivity       (size N_frequencies)
 * @param norm_c_therm_emissivity_epsilon_c Normalized thermal opacity          (size N_frequencies)
 * @param N_frequencies  Number of frequency points
 * @param sigma_c            Output reconstructed scattering opacity  (size N_frequencies)
 * @param epsilon_c          Output reconstructed thermal emissivity  (size N_frequencies)
 * @param c_therm_emissivity_epsilon_c      Output reconstructed thermal opacity     (size N_frequencies)
 * @return 0 on success, -1 on error
 */
int
denormalize_continuum_data(THDF_continuum_t cont, const uint16_t *norm_sigma_c, const uint16_t *norm_epsilon_c,
                           const uint16_t *norm_c_therm_emissivity_epsilon_c, int N_frequencies, double *sigma_c,
                           double *epsilon_c, double *c_therm_emissivity_epsilon_c);

/**
 * @brief Write normalized continuum data for a spatial block to HDF5
 *
 * Writes a [count_x, count_y, count_z, N_frequencies] hyperslab to each
 * 4D uint16 dataset and the normalization struct to the 3D compound dataset.
 *
 * @param file            HDF5 file handle
 * @param N_frequencies   Number of frequency points
 * @param ix, iy, iz      Starting spatial indices
 * @param count_x, count_y, count_z  Block extents in each spatial dimension
 * @param normalized_cont Pointer to normalization scalars for this block
 * @param norm_sigma_c       uint16 scattering opacity    (size N_frequencies)
 * @param norm_epsilon_c     uint16 thermal emissivity    (size N_frequencies)
 * @param norm_c_therm_emissivity_epsilon_c uint16 thermal opacity       (size N_frequencies)
 * @return 0 on success, -1 on error
 */
int
write_normalized_continuum_to_hdf5(hid_t file, int N_frequencies, int ix, int iy, int iz, int count_x, int count_y,
                                   int count_z, const THDF_continuum_t *normalized_cont, uint16_t *norm_sigma_c,
                                   uint16_t *norm_epsilon_c, uint16_t *norm_c_therm_emissivity_epsilon_c);

int
read_continuum_from_hdf5(hid_t     file,                                                               //
                         const int N_frequencies,                                                      //
                         const int ix, const int iy, const int iz,                                     //
                         const int count_x, const int count_y, const int count_z,                      //
                         THDF_continuum_t *cont_to_fill,                                               //
                         double *sigma_c, double *epsilon_c, double *c_therm_emissivity_epsilon_c);  //

/**
 * @brief Standalone continuum demo
 *
 * Opens (or creates) an HDF5 file, writes sample continuum data for
 * several grid points, and closes the file.
 *
 * @param argc Argument count (argv[1] may override the default filename)
 * @param argv Argument vector
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int
continuum_main_demo(int argc, char **argv);

#ifdef __cplusplus
}
#endif

#endif /* HDF_ATMOS_CONTINUUM_H */
