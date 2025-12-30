#ifndef __THDF_OUTPUT_AR_FIELD_H__
#define __THDF_OUTPUT_AR_FIELD_H__

#include <stdio.h>
#include <stdlib.h>
#include "hdf_output_field.h"

#ifdef __cplusplus
extern "C"
{
#endif

	/**
	 * @brief Structure to hold arbitrary field data for a single grid point and beam
	 *
	 */
	typedef struct
	{
		int index_i;	// horizontal x index
		int index_j;	// horizontal y index
		int index_k;	// height index
		int beam_index; // arbitrary beam index

		int N_frequencies;

		/**
		 * @brief Stokes parameters I, Q, U, V for each frequency
		 * @note Stored as 3D ROW major arrays [i, k, frequency] where i and k are spatial indices
		 * and frequency is the frequency index
		 */
		THDF_float_t *stokes_I;	 //
		THDF_float_t *stokes_QI; //
		THDF_float_t *stokes_UI; //
		THDF_float_t *stokes_VI; //

		THDF_n_float_t *norm_multiplier_I;	// TODO: not used yet
		THDF_n_float_t *norm_multiplier_QI; //
		THDF_n_float_t *norm_multiplier_UI; //
		THDF_n_float_t *norm_multiplier_VI; //

	} THDF_ard_field_t;

	/**
	 * @brief HDF5 handler for arbitrary field output datasets.
	 *
	 * This structure manages HDF5 file, group, dataset, and dataspace identifiers
	 * for storing Stokes parameters (I, Q, U, V) for arbitrary beams and frequencies.
	 * It also tracks the grid size and state of the handler.
	 */
	typedef struct
	{
		hid_t file_id;	/**< HDF5 file identifier */
		hid_t group_id; /**< HDF5 group identifier for this field */

		hid_t dataset_id_I; /**< HDF5 dataset identifier for Stokes I */
		hid_t dataset_id_Q; /**< HDF5 dataset identifier for Stokes Q */
		hid_t dataset_id_U; /**< HDF5 dataset identifier for Stokes U */
		hid_t dataset_id_V; /**< HDF5 dataset identifier for Stokes V */

		// hid_t dataset_id_mu;
		// hid_t dataset_id_chi;

		// hid_t dataspace_mu;
		// hid_t dataspace_chi;

		hid_t dataspace_id_I; /**< HDF5 dataspace identifier for Stokes I */
		hid_t dataspace_id_Q; /**< HDF5 dataspace identifier for Stokes Q */
		hid_t dataspace_id_U; /**< HDF5 dataspace identifier for Stokes U */
		hid_t dataspace_id_V; /**< HDF5 dataspace identifier for Stokes V */

		hid_t datatype_id;		 /**< HDF5 datatype identifier for the stored data */
		int	  N_x;				 /**< Number of grid points in x direction */
		int	  N_y;				 /**< Number of grid points in y direction */
		int	  N_arbitrary_beams; /**< Number of arbitrary beams */
		int	  N_frequencies;	 /**< Number of frequencies */
		int	  is_open;			 /**< Handler open state flag (1=open, 0=closed) */
	} THDF_ard_field_handler;

	/**
	 * @brief Structure to hold arbitrary directions (mu, theta)
	 * namely the cosines of the inclination angles and the azimuthal angles.
	 */
	typedef struct
	{
		double *mu;	 /**< Cosines of the inclination angles */
		double *chi; /**< Azimuthal angles */

		int N_directions; /**< Number of arbitrary directions */

	} THDF_ard_directions_t;

	/**
	 * @brief Structure to hold arbitrary frequencies.
	 */
	typedef struct
	{
		double *frequencies;   /**< Array of frequencies */
		int		N_frequencies; /**< Number of frequencies */

	} THDF_ard_frequencies_t;

	/** @brief Build an empty THDF_ard_field_handler structure */
	THDF_ard_field_handler
	THDF_make_empty_ard_field_handler(void);

	/** @brief
	 * Create and initialize an HDF5 handler for arbitrary field output datasets.
	 * This function sets up the necessary HDF5 groups, dataspaces, and datasets
	 * @param file_id HDF5 file identifier
	 * @param mu Cosine of inclination angle (not used in this function)
	 * @param chi Azimuthal angle (not used in this function)
	 * @param N_x Number of grid points in x direction
	 * @param N_y Number of grid points in y direction
	 * @param N_arbitrary_beams Number of arbitrary beams
	 * @param N_frequencies Number of frequencies
	 * @return Pointer to initialized THDF_ard_field_handler structure, or NULL on error
	 */
	THDF_ard_field_handler *
	THDF_create_ard_field_handler_mpi(hid_t	 file_id,			//
									  double mu,				//
									  double chi,				//
									  int	 N_x,				//
									  int	 N_y,				//
									  int	 N_arbitrary_beams, //
									  int	 N_frequencies);		//

	/** @brief
	 * Close and free resources associated with an HDF5 handler for arbitrary field output datasets.
	 * @param output_dset Pointer to THDF_ard_field_handler structure to close
	 */
	void
	THDF_close_ard_field_handler_mpi(THDF_ard_field_handler *output_dset);

	/** @brief
	 * Write arbitrary field data to HDF5 datasets.
	 * @param output_dset Pointer to THDF_ard_field_handler structure
	 * @param output_field Pointer to THDF_ard_field_t structure containing data to write
	 * @param start_i Starting index in x direction
	 * @param start_j Starting index in y direction
	 * @param count_i Number of elements to write in x direction
	 * @param count_j Number of elements to write in y direction
	 * @param count_frequencies Number of frequency elements to write
	 * @return 0 on success, -1 on error
	 */
	int
	THDF_write_ard_field_to_hdf5(THDF_ard_field_handler *output_dset,  //
								 THDF_ard_field_t		*output_field, //
								 const hsize_t			 start_i,	   //
								 const hsize_t			 start_j,	   //
								 const hsize_t			 count_i,	   //
								 const hsize_t			 count_j,	   //
								 const hsize_t			 count_frequencies);	   //

	/** @brief
	 * Write arbitrary directions (mu, chi) to HDF5 file.
	 * @param file_id HDF5 file identifier
	 * @param ard_directions Pointer to THDF_ard_directions_t structure containing directions to write
	 * @return 0 on success, -1 on error
	 */
	int
	THDF_write_ard_directions(hid_t						   file_id,			//
							  const THDF_ard_directions_t *ard_directions); //

	/** @brief
	 * Write arbitrary frequencies to HDF5 file.
	 * @param file_id HDF5 file identifier
	 * @param ard_freqs Pointer to THDF_ard_frequencies_t structure containing frequencies to write
	 * @return 0 on success, -1 on error
	 */
	int
	THDF_write_ard_frequencies(hid_t						 file_id,	 //
							   const THDF_ard_frequencies_t *ard_freqs); //

#ifdef __cplusplus
} // extern "C"
#endif

#endif // __THDF_OUTPUT_AR_FIELD_H__