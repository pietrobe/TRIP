#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf_output_3D_field.h"

/////////////////////////////////////////////////////
// write_stokes_hyperslab
/////////////////////////////////////////////////////
// write_stokes_hyperslab riceve dxpl dall'esterno
static int
write_stokes_hyperslab(hid_t		  dataset_id,  //
					   hid_t		  datatype_id, //
					   hid_t		  memspace,	   //
					   const hsize_t *start,	   //
					   const hsize_t *count,	   //
					   const void	 *data,		   //
					   hid_t		  dxpl,		   // ← ricevuto
					   const char	 *stokes_name)
{
	hid_t dataspace = H5Dget_space(dataset_id);
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, NULL);
	herr_t status = H5Dwrite(dataset_id, datatype_id, memspace, dataspace, dxpl, data); // ← usa quello esterno
	HDF5_CHECK_WRITE_OF(status, stokes_name, memspace, -1);
	H5Sclose(dataspace);
	return (int)status;
}
/////////////////////////////////////////////////////
// get_input_index
/////////////////////////////////////////////////////
static inline int
get_input_index(const int i,					 //
				const int j,					 //
				const int k,					 //
				const int incl,					 //
				const int az,					 //
				const int f,					 //
				const int stokes_offset,		 //
				const int N0_in_local_x,		 //
				const int N0_in_local_y,		 //
				const int N0_in_local_z,		 //
				const int stride_in_x,			 //
				const int stride_in_y,			 //
				const int stride_in_z,			 //
				const int stride_in_incl,		 //
				const int stride_in_azimuth,	 //
				const int stride_in_frequencies, //
				const int stride_in_stokes)
{ //

	return (i + N0_in_local_x) * stride_in_x + (j + N0_in_local_y) * stride_in_y + (k + N0_in_local_z) * stride_in_z
		   + incl * stride_in_incl + az * stride_in_azimuth + f * stride_in_frequencies
		   + stokes_offset * stride_in_stokes;
}

/*
 * THDF_create_3D_field_handler_mpi
 *
 * Crea il file HDF5 e la struttura dei dataset da rank 0 in modo seriale.
 * Ottimizzata per minimizzare il tempo di creazione su filesystem paralleli:
 *
 *   - meta_block_size grande: aggrega tutti i metadata write in un
 *     unico blocco contiguo invece di tante piccole scritture sparse
 *   - libver LATEST: abilita v2 B-tree e object header v2, più compatti
 *   - file_space PAGE: layout a pagine, metadata su un numero minimo
 *     di OST stripe unit
 *   - small_data_block_size = 0: nessuna aggregazione dati (non serve
 *     per una funzione che scrive solo struttura)
 *   - tutti i dataset creati in sequenza con lo stesso fspace/dcpl
 *     riusato: un solo H5Screate_simple per dimensionalità
 *   - H5Fflush(GLOBAL) esplicito prima di H5Fclose: garantisce che
 *     il filesystem abbia ricevuto tutto prima della MPI_Barrier
 *     del chiamante
 *
 * Ritorna 0 in caso di successo, -1 in caso di errore.
 *
 * Uso tipico:
 *   if (rank == 0)
 *       if (THDF_create_3D_field_handler_mpi(path, ...) < 0)
 *           MPI_Abort(comm, 1);
 *   MPI_Barrier(comm);
 *   // tutti i writer aprono con THDF_open_3D_field_handler_mpi
 */

#define HDF_CLOSE_IF(id, fn) \
	do                       \
	{                        \
		if ((id) >= 0)       \
		{                    \
			fn(id);          \
			(id) = -1;       \
		}                    \
	} while (0)

int
THDF_create_3D_field_handler_mpi(const char *path, const bool normalized_output, const int N_x, const int N_y,
								 const int N_z, const int N_incl, const int N_azimuth, const int N_frequencies)
{
	hid_t fapl	 = -1;
	hid_t file	 = -1;
	hid_t group	 = -1;
	hid_t dcpl	 = -1;
	hid_t dcpl_n = -1;
	hid_t fspace = -1;
	hid_t nspace = -1;
	hid_t f32	 = -1;
	hid_t f64	 = -1;
	hid_t dset	 = -1;
	int	  rc	 = 0;

	/* ----------------------------------------------------------------
	 * FAPL minimale — solo meta_block_size.
	 * H5Pset_libver_bounds con H5F_LIBVER_V18 è invalido in HDF5 2.0:
	 * la funzione fallisce silenziosamente corrompendo il fapl.
	 * Rimosso. H5F_LIBVER_LATEST come singolo bound è applicabile
	 * solo se entrambi low/high sono LATEST — non misto con V18.
	 * ---------------------------------------------------------------- */
	fapl = H5Pcreate(H5P_FILE_ACCESS);
	if (fapl < 0)
	{
		fprintf(stderr, "[rank 0] H5Pcreate(FILE_ACCESS) fallito\n");
		return -1;
	}

	/*
	 * Aggrega tutti i metadata write in un unico blocco contiguo.
	 * Scaricato al flush esplicito prima del close.
	 */
	if (H5Pset_meta_block_size(fapl, 4 * 1024 * 1024) < 0)
	{
		fprintf(stderr, "[rank 0] H5Pset_meta_block_size fallito\n");
		rc = -1;
		goto cleanup;
	}

	/* ---- Crea il file ---- */
	file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
	if (file < 0)
	{
		fprintf(stderr, "[rank 0] Impossibile creare il file: %s\n", path);
		rc = -1;
		goto cleanup;
	}

	/* fapl non serve più */
	HDF_CLOSE_IF(fapl, H5Pclose);

	/* ---- Gruppo ---- */
	group = H5Gcreate2(file, "/radiation_field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (group < 0)
	{
		fprintf(stderr, "[rank 0] Errore creando /radiation_field\n");
		rc = -1;
		goto cleanup;
	}

	/* ---- Tipo f32 ---- */
	f32 = H5Tcopy(THDF_get_hdf_float32_datatype());
	if (f32 < 0)
	{
		rc = -1;
		goto cleanup;
	}

	/* ---- dcpl 6D ---- */
	dcpl = H5Pcreate(H5P_DATASET_CREATE);
	if (dcpl < 0)
	{
		rc = -1;
		goto cleanup;
	}
	{
		hsize_t chunk6[6] = {1, (hsize_t)N_y, (hsize_t)N_z, (hsize_t)N_incl, (hsize_t)N_azimuth, (hsize_t)N_frequencies};
		if (H5Pset_chunk(dcpl, 6, chunk6) < 0)
		{
			rc = -1;
			goto cleanup;
		}
		if (H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY) < 0)
		{
			rc = -1;
			goto cleanup;
		}
	}

	/* ---- Dataset 6D ---- */
	{
		hsize_t dims6[6] = {(hsize_t)N_x,	 (hsize_t)N_y,		 (hsize_t)N_z,
							(hsize_t)N_incl, (hsize_t)N_azimuth, (hsize_t)N_frequencies};
		fspace			 = H5Screate_simple(6, dims6, NULL);
		if (fspace < 0)
		{
			rc = -1;
			goto cleanup;
		}

		const char *names6[4] = {"I", "QI_pc", "UI_pc", "VI_pc"};
		for (int i = 0; i < 4; i++)
		{
			dset = H5Dcreate2(group, names6[i], f32, fspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);
			if (dset < 0)
			{
				fprintf(stderr, "[rank 0] Errore creando dataset %s\n", names6[i]);
				rc = -1;
				goto cleanup;
			}
			HDF_CLOSE_IF(dset, H5Dclose);
		}
		HDF_CLOSE_IF(fspace, H5Sclose);
	}
	HDF_CLOSE_IF(dcpl, H5Pclose);

	if (!normalized_output) goto sync;

	/* ---- Tipo f64 ---- */
	f64 = H5Tcopy(THDF_get_hdf_float_datatype());
	if (f64 < 0)
	{
		rc = -1;
		goto cleanup;
	}

	/* ---- dcpl 5D ---- */
	dcpl_n = H5Pcreate(H5P_DATASET_CREATE);
	if (dcpl_n < 0)
	{
		rc = -1;
		goto cleanup;
	}
	{
		hsize_t chunk5[5] = {1, (hsize_t)N_y, (hsize_t)N_z, (hsize_t)N_incl, (hsize_t)N_azimuth};
		if (H5Pset_chunk(dcpl_n, 5, chunk5) < 0)
		{
			rc = -1;
			goto cleanup;
		}
		if (H5Pset_alloc_time(dcpl_n, H5D_ALLOC_TIME_EARLY) < 0)
		{
			rc = -1;
			goto cleanup;
		}
	}

	/* ---- Dataset 5D ---- */
	{
		hsize_t dims5[5] = {(hsize_t)N_x, (hsize_t)N_y, (hsize_t)N_z, (hsize_t)N_incl, (hsize_t)N_azimuth};
		nspace			 = H5Screate_simple(5, dims5, NULL);
		if (nspace < 0)
		{
			rc = -1;
			goto cleanup;
		}

		const char *names5[4] = {"norm_multiplier_I", "norm_multiplier_QI_pc", "norm_multiplier_UI_pc",
								 "norm_multiplier_VI_pc"};
		for (int i = 0; i < 4; i++)
		{
			dset = H5Dcreate2(group, names5[i], f64, nspace, H5P_DEFAULT, dcpl_n, H5P_DEFAULT);
			if (dset < 0)
			{
				fprintf(stderr, "[rank 0] Errore creando dataset %s\n", names5[i]);
				rc = -1;
				goto cleanup;
			}
			HDF_CLOSE_IF(dset, H5Dclose);
		}
		HDF_CLOSE_IF(nspace, H5Sclose);
	}
	HDF_CLOSE_IF(dcpl_n, H5Pclose);

sync:
	/*
	 * Flush esplicito: svuota il meta_block_size buffer verso il
	 * filesystem prima che gli altri rank aprano il file.
	 * H5Fclose da solo non è sufficiente su filesystem con
	 * writeback caching.
	 */
	if (H5Fflush(file, H5F_SCOPE_GLOBAL) < 0)
	{
		fprintf(stderr, "[rank 0] H5Fflush fallito — il file potrebbe essere incompleto\n");
		rc = -1;
	}

cleanup:
	HDF_CLOSE_IF(dset, H5Dclose);
	HDF_CLOSE_IF(fspace, H5Sclose);
	HDF_CLOSE_IF(nspace, H5Sclose);
	HDF_CLOSE_IF(dcpl, H5Pclose);
	HDF_CLOSE_IF(dcpl_n, H5Pclose);
	HDF_CLOSE_IF(f32, H5Tclose);
	HDF_CLOSE_IF(f64, H5Tclose);
	HDF_CLOSE_IF(group, H5Gclose);
	HDF_CLOSE_IF(fapl, H5Pclose);

	if (file >= 0) H5Fclose(file);

	return rc;
}

/*
 * THDF_open_3D_field_handler_mpi
 *
 * Apre collettivamente una struttura HDF5 già creata da rank 0.
 * Stessa firma di THDF_create_3D_field_handler_mpi, solo H5Dopen2
 * al posto di H5Dcreate2 — nessuna allocazione su disco.
 *
 * Precondizioni (chiamante):
 *   - rank 0 ha già creato il file con THDF_create_3D_field_handler_mpi
 *     (o THDF_create_3D_field_structure) e lo ha chiuso
 *   - MPI_Barrier eseguita prima di chiamare questa funzione
 *   - `file` aperto con H5Pset_fapl_mpio su writer_comm
 *   - H5Pset_all_coll_metadata_ops + H5Pset_coll_metadata_write attivi
 */
THDF_3D_field_handler_t *
THDF_open_3D_field_handler_mpi(hid_t file, const bool normalized_output, const int N_x, const int N_y, const int N_z,
							   const int N_incl, const int N_azimuth, const int N_frequencies)
{
	/* ---- Allocazione e inizializzazione ---- */
	THDF_3D_field_handler_t *h = malloc(sizeof(*h));
	if (!h)
	{
		fprintf(stderr, "Error allocating memory for output field dataset struct\n");
		return NULL;
	}

	h->file_id						= file;
	h->group_id						= -1;
	h->dataset_id_I					= -1;
	h->dataset_id_Q					= -1;
	h->dataset_id_U					= -1;
	h->dataset_id_V					= -1;
	h->dataspace_id_I				= -1;
	h->dataspace_id_Q				= -1;
	h->dataspace_id_U				= -1;
	h->dataspace_id_V				= -1;
	h->dataset_norm_multiplier_I	= -1;
	h->dataset_norm_multiplier_QI	= -1;
	h->dataset_norm_multiplier_UI	= -1;
	h->dataset_norm_multiplier_VI	= -1;
	h->dataspace_norm_multiplier_I	= -1;
	h->dataspace_norm_multiplier_QI = -1;
	h->dataspace_norm_multiplier_UI = -1;
	h->dataspace_norm_multiplier_VI = -1;
	h->datatype_id					= -1;
	h->datatype_f32_id				= -1;
	h->is_open						= false;

	h->N_x				 = N_x;
	h->N_y				 = N_y;
	h->N_z				 = N_z;
	h->N_incl			 = N_incl;
	h->N_azimuth		 = N_azimuth;
	h->N_frequencies	 = N_frequencies;
	h->normalized_output = normalized_output;

	/* ---- Gruppo ---- */
	/*
	 * H5Gopen2 è collettivo con H5Pset_all_coll_metadata_ops attivo.
	 * Nessun H5E_BEGIN_TRY: il gruppo deve esistere per definizione,
	 * se non esiste è un errore del chiamante.
	 */
	h->group_id = H5Gopen2(file, "/radiation_field", H5P_DEFAULT);
	if (h->group_id < 0)
	{
		fprintf(stderr, "Error opening /radiation_field group\n");
		goto fail;
	}

	/* ---- Tipi ---- */
	h->datatype_f32_id = H5Tcopy(THDF_get_hdf_float32_datatype());
	if (h->datatype_f32_id < 0)
	{
		goto fail;
	}

	/* ---- Dataset 6D ---- */
	{
		const char *names6[4]  = {"I", "QI_pc", "UI_pc", "VI_pc"};
		hid_t	   *targets[4] = {&h->dataset_id_I, &h->dataset_id_Q, &h->dataset_id_U, &h->dataset_id_V};

		for (int i = 0; i < 4; i++)
		{
			*targets[i] = H5Dopen2(h->group_id, names6[i], H5P_DEFAULT);
			if (*targets[i] < 0)
			{
				fprintf(stderr, "Error opening dataset %s\n", names6[i]);
				goto fail;
			}
		}
	}

	if (!normalized_output)
	{
		h->is_open = true;
		return h;
	}

	/* ---- Tipo f64 per i moltiplicatori ---- */
	h->datatype_id = H5Tcopy(THDF_get_hdf_float_datatype());
	if (h->datatype_id < 0)
	{
		goto fail;
	}

	/* ---- Dataset 5D (moltiplicatori di normalizzazione) ---- */
	{
		const char *names5[4]		= {"norm_multiplier_I", "norm_multiplier_QI_pc", "norm_multiplier_UI_pc",
									   "norm_multiplier_VI_pc"};
		hid_t	   *norm_targets[4] = {&h->dataset_norm_multiplier_I, &h->dataset_norm_multiplier_QI,
									   &h->dataset_norm_multiplier_UI, &h->dataset_norm_multiplier_VI};

		for (int i = 0; i < 4; i++)
		{
			*norm_targets[i] = H5Dopen2(h->group_id, names5[i], H5P_DEFAULT);
			if (*norm_targets[i] < 0)
			{
				fprintf(stderr, "Error opening dataset %s\n", names5[i]);
				goto fail;
			}
		}
	}

	h->is_open = true;
	return h;

fail:
	THDF_close_3D_field_handler_mpi(h);
	return NULL;
}

/////////////////////////////////////////////////////
// THDF_close_3D_field_handler_mpi
/////////////////////////////////////////////////////
void
THDF_close_3D_field_handler_mpi(THDF_3D_field_handler_t *output_dset)
{
	if (!output_dset)
	{
		return;
	}

	if (output_dset->is_open)
	{
		if (output_dset->dataset_id_I >= 0)
		{
			H5Dclose(output_dset->dataset_id_I);
		}
		if (output_dset->dataset_id_Q >= 0)
		{
			H5Dclose(output_dset->dataset_id_Q);
		}
		if (output_dset->dataset_id_U >= 0)
		{
			H5Dclose(output_dset->dataset_id_U);
		}
		if (output_dset->dataset_id_V >= 0)
		{
			H5Dclose(output_dset->dataset_id_V);
		}
		if (output_dset->dataspace_id_I >= 0)
		{
			H5Sclose(output_dset->dataspace_id_I);
		}
		if (output_dset->dataspace_id_Q >= 0)
		{
			H5Sclose(output_dset->dataspace_id_Q);
		}
		if (output_dset->dataspace_id_U >= 0)
		{
			H5Sclose(output_dset->dataspace_id_U);
		}
		if (output_dset->dataspace_id_V >= 0)
		{
			H5Sclose(output_dset->dataspace_id_V);
		}

		if (output_dset->dataset_norm_multiplier_I >= 0)
		{
			H5Dclose(output_dset->dataset_norm_multiplier_I);
		}
		if (output_dset->dataset_norm_multiplier_QI >= 0)
		{
			H5Dclose(output_dset->dataset_norm_multiplier_QI);
		}
		if (output_dset->dataset_norm_multiplier_UI >= 0)
		{
			H5Dclose(output_dset->dataset_norm_multiplier_UI);
		}
		if (output_dset->dataset_norm_multiplier_VI >= 0)
		{
			H5Dclose(output_dset->dataset_norm_multiplier_VI);
		}
		if (output_dset->dataspace_norm_multiplier_I >= 0)
		{
			H5Sclose(output_dset->dataspace_norm_multiplier_I);
		}
		if (output_dset->dataspace_norm_multiplier_QI >= 0)
		{
			H5Sclose(output_dset->dataspace_norm_multiplier_QI);
		}
		if (output_dset->dataspace_norm_multiplier_UI >= 0)
		{
			H5Sclose(output_dset->dataspace_norm_multiplier_UI);
		}
		if (output_dset->dataspace_norm_multiplier_VI >= 0)
		{
			H5Sclose(output_dset->dataspace_norm_multiplier_VI);
		}

		if (output_dset->datatype_id >= 0)
		{
			H5Tclose(output_dset->datatype_id);
		}
		if (output_dset->datatype_f32_id >= 0)
		{
			H5Tclose(output_dset->datatype_f32_id);
		}

		if (output_dset->group_id >= 0)
		{
			H5Gclose(output_dset->group_id);
		}
		output_dset->is_open = false;
	}

	free(output_dset);
}

/////////////////////////////////////////////////////
// THDF_write_3D_field_dataset_to_hdf5
/////////////////////////////////////////////////////
int
THDF_write_3D_field_dataset_to_hdf5(THDF_3D_field_handler_t *output_dset, THDF_3D_field_t *output_field, hsize_t start_i,
									hsize_t start_j, hsize_t start_k, hsize_t start_incl, hsize_t start_azimuth,
									hsize_t count_i, hsize_t count_j, hsize_t count_k, hsize_t count_incl,
									hsize_t count_azimuth, hsize_t count_frequencies)
{
	if (!output_dset || !output_field || !output_dset->is_open)
	{
		fprintf(stderr, "Invalid handler\n");
		return -1;
	}

	/* ---- DXPL collettivo ---- */
	hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
	if (dxpl < 0) return -1;
	H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

	hsize_t start[6] = {start_i, start_j, start_k, start_incl, start_azimuth, 0};
	hsize_t count[6] = {count_i, count_j, count_k, count_incl, count_azimuth, count_frequencies};

	hid_t memspace = H5Screate_simple(6, count, NULL);

	write_stokes_hyperslab(output_dset->dataset_id_I, output_dset->datatype_f32_id, memspace, start, count,
						   output_field->stokes_I, dxpl, "Stokes I");
	write_stokes_hyperslab(output_dset->dataset_id_Q, output_dset->datatype_f32_id, memspace, start, count,
						   output_field->stokes_QI, dxpl, "Stokes QI");
	write_stokes_hyperslab(output_dset->dataset_id_U, output_dset->datatype_f32_id, memspace, start, count,
						   output_field->stokes_UI, dxpl, "Stokes UI");
	write_stokes_hyperslab(output_dset->dataset_id_V, output_dset->datatype_f32_id, memspace, start, count,
						   output_field->stokes_VI, dxpl, "Stokes VI");
	H5Sclose(memspace);

/* ---- verifica che la collective I/O non sia degradata ---- */
#ifndef NDEBUG
	{
		uint32_t local_cause, global_cause;
		H5Pget_mpio_no_collective_cause(dxpl, &local_cause, &global_cause);
		if (global_cause != 0)
		{
			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			if (rank == 0) fprintf(stderr, "WARN: collective I/O fallback cause=0x%x\n", global_cause);
		}
	}
#endif

	if (output_dset->normalized_output)
	{
		hsize_t norm_start[5] = {start_i, start_j, start_k, start_incl, start_azimuth};
		hsize_t norm_count[5] = {count_i, count_j, count_k, count_incl, count_azimuth};
		hid_t	norm_memspace = H5Screate_simple(5, norm_count, NULL);

		write_stokes_hyperslab(output_dset->dataset_norm_multiplier_I, output_dset->datatype_id, norm_memspace,
							   norm_start, norm_count, output_field->norm_multiplier_I, dxpl, "norm I");
		write_stokes_hyperslab(output_dset->dataset_norm_multiplier_QI, output_dset->datatype_id, norm_memspace,
							   norm_start, norm_count, output_field->norm_multiplier_QI, dxpl, "norm QI");
		write_stokes_hyperslab(output_dset->dataset_norm_multiplier_UI, output_dset->datatype_id, norm_memspace,
							   norm_start, norm_count, output_field->norm_multiplier_UI, dxpl, "norm UI");
		write_stokes_hyperslab(output_dset->dataset_norm_multiplier_VI, output_dset->datatype_id, norm_memspace,
							   norm_start, norm_count, output_field->norm_multiplier_VI, dxpl, "norm VI");
		H5Sclose(norm_memspace);
	}

	H5Pclose(dxpl);
	return 0;
}

/////////////////////////////////////////////////////
// THDF_copy_3D_block_field
/////////////////////////////////////////////////////
static inline double
read_stokes_input_value(const void *stokes_IQUI, const bool input_is_float32, const int idx)
{
	if (input_is_float32)
	{
		return (double)((const THDF_float32_t *)stokes_IQUI)[idx];
	}

	return (double)((const THDF_float_t *)stokes_IQUI)[idx];
}

static void
copy_3D_block_field_impl(THDF_3D_field_t *field,			 //
						 const bool		  normalized_output, //
						 const void		 *stokes_IQUI,		 //
						 const bool		  input_is_float32,	 //
						 THDF_float32_t	 *stokes_out_I,		 //
						 THDF_float32_t *stokes_out_QI, THDF_float32_t *stokes_out_UI, THDF_float32_t *stokes_out_VI,
						 THDF_float_t *norm_multiplier_I, THDF_float_t *norm_multiplier_QI,
						 THDF_float_t *norm_multiplier_UI, THDF_float_t *norm_multiplier_VI, const int N0_in_local_x,
						 const int N0_in_local_y, const int N0_in_local_z, const int N_local_x, const int N_local_y,
						 const int N_local_z, const int N_incl, const int N_azimuth, const int N_frequencies,
						 const int stride_in_x, const int stride_in_y, const int stride_in_z, const int stride_in_incl,
						 const int stride_in_azimuth, const int stride_in_frequencies, const int stride_in_stokes)
{
	if (normalized_output)
	{
		field->norm_multiplier_I  = norm_multiplier_I;
		field->norm_multiplier_QI = norm_multiplier_QI;
		field->norm_multiplier_UI = norm_multiplier_UI;
		field->norm_multiplier_VI = norm_multiplier_VI;
	}
	else
	{
		field->norm_multiplier_I  = NULL;
		field->norm_multiplier_QI = NULL;
		field->norm_multiplier_UI = NULL;
		field->norm_multiplier_VI = NULL;
	}

	field->stokes_I	 = stokes_out_I;
	field->stokes_QI = stokes_out_QI;
	field->stokes_UI = stokes_out_UI;
	field->stokes_VI = stokes_out_VI;

	for (int i = 0; i < N_local_x; i++)
	{
		for (int j = 0; j < N_local_y; j++)
		{
			for (int k = 0; k < N_local_z; k++)
			{
				for (int incl = 0; incl < N_incl; incl++)
				{
					for (int az = 0; az < N_azimuth; az++)
					{
						int norm_idx = ((((i * N_local_y + j) * N_local_z + k) * N_incl + incl) * N_azimuth) + az;

						double abs_max_I  = 0.0;
						double abs_max_QI = 0.0;
						double abs_max_UI = 0.0;
						double abs_max_VI = 0.0;

						if (normalized_output)
						{
							for (int f = 0; f < N_frequencies; f++)
							{
								const int in_index_I = get_input_index(
									i, j, k, incl, az, f, 0, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
									stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
									stride_in_stokes); // Stokes I is at offset 0 in the stokes dimension

								const int in_index_QI = get_input_index(
									i, j, k, incl, az, f, 1, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
									stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
									stride_in_stokes); // Stokes Q/I is at offset 1 in the stokes dimension

								const int in_index_UI = get_input_index(
									i, j, k, incl, az, f, 2, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
									stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
									stride_in_stokes); // Stokes U/I is at offset 2 in the stokes dimension

								const int in_index_VI = get_input_index(
									i, j, k, incl, az, f, 3, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
									stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
									stride_in_stokes); // Stokes V/I is at offset 3 in the stokes dimension

								const double val_I	= read_stokes_input_value(stokes_IQUI, input_is_float32, in_index_I);
								const double val_QI = read_stokes_input_value(stokes_IQUI, input_is_float32, in_index_QI);
								const double val_UI = read_stokes_input_value(stokes_IQUI, input_is_float32, in_index_UI);
								const double val_VI = read_stokes_input_value(stokes_IQUI, input_is_float32, in_index_VI);

								if (fabs(val_I) > abs_max_I)
								{
									abs_max_I = fabs(val_I);
								}

								if (fabs(val_QI) > abs_max_QI)
								{
									abs_max_QI = fabs(val_QI);
								}

								if (fabs(val_UI) > abs_max_UI)
								{
									abs_max_UI = fabs(val_UI);
								}

								if (fabs(val_VI) > abs_max_VI)
								{
									abs_max_VI = fabs(val_VI);
								}
							}
						}

						const double norm_multiplier_I_	 = (abs_max_I > 0.0) ? abs_max_I : 1.0;
						const double norm_multiplier_QI_ = (abs_max_QI > 0.0) ? abs_max_QI : 1.0;
						const double norm_multiplier_UI_ = (abs_max_UI > 0.0) ? abs_max_UI : 1.0;
						const double norm_multiplier_VI_ = (abs_max_VI > 0.0) ? abs_max_VI : 1.0;

						if (normalized_output)
						{
							norm_multiplier_I[norm_idx]	 = norm_multiplier_I_;
							norm_multiplier_QI[norm_idx] = norm_multiplier_QI_;
							norm_multiplier_UI[norm_idx] = norm_multiplier_UI_;
							norm_multiplier_VI[norm_idx] = norm_multiplier_VI_;
						}

						for (int f = 0; f < N_frequencies; f++)
						{
							int idx = ((((((i)*N_local_y + (j)) * N_local_z + (k)) * //
											 N_incl
										 + incl)
											* N_azimuth
										+ az)
									   * N_frequencies)
									  + f;

							const int in_index_I = get_input_index(
								i, j, k, incl, az, f, 0, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
								stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
								stride_in_stokes); // Stokes I is at offset 0 in the stokes dimension

							const int in_index_QI = get_input_index(
								i, j, k, incl, az, f, 1, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
								stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
								stride_in_stokes); // Stokes Q/I is at offset 1 in the stokes dimension

							const int in_index_UI = get_input_index(
								i, j, k, incl, az, f, 2, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
								stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
								stride_in_stokes); // Stokes U/I is at offset 2 in the stokes dimension

							const int in_index_VI = get_input_index(
								i, j, k, incl, az, f, 3, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
								stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
								stride_in_stokes); // Stokes V/I is at offset 3 in the stokes dimension

							const double II = read_stokes_input_value(stokes_IQUI, input_is_float32, in_index_I);
							const double QI =
								read_stokes_input_value(stokes_IQUI, input_is_float32, in_index_QI) / II * 100.0;
							const double UI =
								read_stokes_input_value(stokes_IQUI, input_is_float32, in_index_UI) / II * 100.0;
							const double VI =
								read_stokes_input_value(stokes_IQUI, input_is_float32, in_index_VI) / II * 100.0;

							field->stokes_I[idx]  = (THDF_float32_t)II;
							field->stokes_QI[idx] = (THDF_float32_t)QI;
							field->stokes_UI[idx] = (THDF_float32_t)UI;
							field->stokes_VI[idx] = (THDF_float32_t)VI;
						}
					}
				}
			}
		}
	}
}

void
THDF_copy_3D_block_field(THDF_3D_field_t *field,			 //
						 const bool		  normalized_output, //
						 THDF_float_t	 *stokes_IQUI,		 // Input data in float 64 bits
						 THDF_float32_t	 *stokes_out_I,		 // Output normalized buffer data in float 32 bits
						 THDF_float32_t	 *stokes_out_QI,	 // This is a preallocated buffer that
						 THDF_float32_t	 *stokes_out_UI,	 // will hold the normalized Stokes Q/I values
						 THDF_float32_t	 *stokes_out_VI,	 //
						 THDF_float_t *norm_multiplier_I,  // Preallocated buffers that will hold the output normalization
						 THDF_float_t *norm_multiplier_QI, // multipliers for each Stokes parameter
						 THDF_float_t *norm_multiplier_UI, // Set to NULL if normalized_output is false
						 THDF_float_t *norm_multiplier_VI, //
						 const int	   N0_in_local_x,	   //
						 const int	   N0_in_local_y,	   //
						 const int	   N0_in_local_z,	   //
						 const int	   N_local_x,		   //
						 const int	   N_local_y,		   //
						 const int	   N_local_z,		   //
						 const int	   N_incl,			   //
						 const int	   N_azimuth,		   //
						 const int	   N_frequencies,	   //
						 const int stride_in_x, // Strides for the input data layout, allowing for flexible arrangements
						 const int stride_in_y, // The function will compute the correct input index using
						 const int stride_in_z, // these strides to access the data in the input arrays
						 const int stride_in_incl,		  //
						 const int stride_in_azimuth,	  //
						 const int stride_in_frequencies, //
						 const int stride_in_stokes)
{
	copy_3D_block_field_impl(field, normalized_output, stokes_IQUI, false, stokes_out_I, stokes_out_QI, stokes_out_UI,
							 stokes_out_VI, norm_multiplier_I, norm_multiplier_QI, norm_multiplier_UI, norm_multiplier_VI,
							 N0_in_local_x, N0_in_local_y, N0_in_local_z, N_local_x, N_local_y, N_local_z, N_incl,
							 N_azimuth, N_frequencies, stride_in_x, stride_in_y, stride_in_z, stride_in_incl,
							 stride_in_azimuth, stride_in_frequencies, stride_in_stokes);
}

void
THDF_copy_3D_block_field_f32(
	THDF_3D_field_t *field,					//
	const bool		 normalized_output,		//
	THDF_float32_t	*stokes_IQUI,			// Input data in float 32 bits
	THDF_float32_t	*stokes_out_I,			// Output normalized buffer data in float 32 bits
	THDF_float32_t	*stokes_out_QI,			// This is a preallocated buffer that
	THDF_float32_t	*stokes_out_UI,			// will hold the normalized Stokes Q/I values
	THDF_float32_t	*stokes_out_VI,			//
	THDF_float_t	*norm_multiplier_I,		// Preallocated buffers that will hold the output normalization
	THDF_float_t	*norm_multiplier_QI,	// multipliers for each Stokes parameter
	THDF_float_t	*norm_multiplier_UI,	// Set to NULL if normalized_output is false
	THDF_float_t	*norm_multiplier_VI,	//
	const int		 N0_in_local_x,			//
	const int		 N0_in_local_y,			//
	const int		 N0_in_local_z,			//
	const int		 N_local_x,				//
	const int		 N_local_y,				//
	const int		 N_local_z,				//
	const int		 N_incl,				//
	const int		 N_azimuth,				//
	const int		 N_frequencies,			//
	const int		 stride_in_x,			// Strides for the input data layout, allowing for flexible arrangements
	const int		 stride_in_y,			// The function will compute the correct input index using
	const int		 stride_in_z,			// these strides to access the data in the input arrays
	const int		 stride_in_incl,		//
	const int		 stride_in_azimuth,		//
	const int		 stride_in_frequencies, //
	const int		 stride_in_stokes)
{
	copy_3D_block_field_impl(field, normalized_output, stokes_IQUI, true, stokes_out_I, stokes_out_QI, stokes_out_UI,
							 stokes_out_VI, norm_multiplier_I, norm_multiplier_QI, norm_multiplier_UI, norm_multiplier_VI,
							 N0_in_local_x, N0_in_local_y, N0_in_local_z, N_local_x, N_local_y, N_local_z, N_incl,
							 N_azimuth, N_frequencies, stride_in_x, stride_in_y, stride_in_z, stride_in_incl,
							 stride_in_azimuth, stride_in_frequencies, stride_in_stokes);
}

//////////////////////////////////////////////////////
// write_3d_field_block_to_hdf5
//////////////////////////////////////////////////////
static void
write_3d_field_block_to_hdf5(hid_t file_id, MPI_Comm comm, bool normalized_output, int N_x, int N_y, int N_z, int N_incl,
							 int N_azimuth, int N_frequencies, int N_local_x, int N_local_y, int N_local_z,
							 int local_start_x, int local_start_y, int local_start_z, THDF_3D_field_t *output_field)
{
	// hid_t file_id = THDF_open_file_MPI(filename, comm);

	// if (file_id < 0) {
	//   MPI_Abort(comm, -1);
	// }

	int comm_rank;
	MPI_Comm_rank(comm, &comm_rank);
	const double write_start = MPI_Wtime();

	THDF_3D_field_handler_t *field_handler = THDF_open_3D_field_handler_mpi(file_id,						   //
																			normalized_output,				   //
																			N_x, N_y, N_z,					   //
																			N_incl, N_azimuth, N_frequencies); //

	THDF_write_3D_field_dataset_to_hdf5(field_handler, output_field, local_start_x, local_start_y, local_start_z, 0, 0,
										N_local_x, N_local_y, N_local_z, N_incl, N_azimuth, N_frequencies);

	const double write_elapsed = MPI_Wtime() - write_start;
	if (comm_rank == 0)
	{
		printf("MPI write elapsed time: %.6f s\n", write_elapsed);
	}

	THDF_close_3D_field_handler_mpi(field_handler);
	// H5Fclose(file_id);
}

//////////////////////////////////////////////////////
// write_3d_field_block_mpi
//////////////////////////////////////////////////////
void
write_3d_field_block_mpi(hid_t			 file_id,			 //
						 MPI_Comm		 comm,				 //
						 bool			 normalized_output,	 //
						 int			 N_x,				 //
						 int			 N_y,				 //
						 int			 N_z,				 //
						 int			 N_incl,			 //
						 int			 N_azimuth,			 //
						 int			 N_frequencies,		 //
						 int			 N_stokes,			 //
						 int			 N_local_x,			 //
						 int			 N_local_y,			 //
						 int			 N_local_z,			 //
						 int			 local_start_x,		 //
						 int			 local_start_y,		 //
						 int			 local_start_z,		 //
						 double			*stokes_IQUI,		 //
						 THDF_float32_t *stokes_I,			 //
						 THDF_float32_t *stokes_QI,			 //
						 THDF_float32_t *stokes_UI,			 //
						 THDF_float32_t *stokes_VI,			 //
						 THDF_float_t	*norm_multiplier_I,	 //
						 THDF_float_t	*norm_multiplier_QI, //
						 THDF_float_t	*norm_multiplier_UI, //
						 THDF_float_t	*norm_multiplier_VI, //
						 int			 stride_x,			 //
						 int			 stride_y,			 //
						 int			 stride_z,			 //
						 int			 stride_incl,		 //
						 int			 stride_azimuth,	 //
						 int			 stride_frequencies, //
						 int			 stride_stokes)
{ //
  //
	(void)N_stokes;

	THDF_3D_field_t output_field;

	THDF_copy_3D_block_field(&output_field, normalized_output, stokes_IQUI, stokes_I, stokes_QI, stokes_UI, stokes_VI,
							 norm_multiplier_I, norm_multiplier_QI, norm_multiplier_UI, norm_multiplier_VI, 0, 0, 0,
							 N_local_x, N_local_y, N_local_z, N_incl, N_azimuth, N_frequencies, stride_x, stride_y,
							 stride_z, stride_incl, stride_azimuth, stride_frequencies, stride_stokes);

	write_3d_field_block_to_hdf5(file_id, comm, normalized_output, N_x, N_y, N_z, N_incl, N_azimuth, N_frequencies,
								 N_local_x, N_local_y, N_local_z, local_start_x, local_start_y, local_start_z,
								 &output_field);
}

void
write_3d_field_block_mpi_f32(hid_t			 file_id,			 //
							 MPI_Comm		 comm,				 //
							 bool			 normalized_output,	 //
							 int			 N_x,				 //
							 int			 N_y,				 //
							 int			 N_z,				 //
							 int			 N_incl,			 //
							 int			 N_azimuth,			 //
							 int			 N_frequencies,		 //
							 int			 N_stokes,			 //
							 int			 N_local_x,			 //
							 int			 N_local_y,			 //
							 int			 N_local_z,			 //
							 int			 local_start_x,		 //
							 int			 local_start_y,		 //
							 int			 local_start_z,		 //
							 THDF_float32_t *stokes_IQUI,		 //
							 THDF_float32_t *stokes_I,			 //
							 THDF_float32_t *stokes_QI,			 //
							 THDF_float32_t *stokes_UI,			 //
							 THDF_float32_t *stokes_VI,			 //
							 THDF_float_t	*norm_multiplier_I,	 //
							 THDF_float_t	*norm_multiplier_QI, //
							 THDF_float_t	*norm_multiplier_UI, //
							 THDF_float_t	*norm_multiplier_VI, //
							 int			 stride_x,			 //
							 int			 stride_y,			 //
							 int			 stride_z,			 //
							 int			 stride_incl,		 //
							 int			 stride_azimuth,	 //
							 int			 stride_frequencies, //
							 int			 stride_stokes)
{ //
	(void)N_stokes;

	THDF_3D_field_t output_field;

	THDF_copy_3D_block_field_f32(&output_field, normalized_output, stokes_IQUI, stokes_I, stokes_QI, stokes_UI, stokes_VI,
								 norm_multiplier_I, norm_multiplier_QI, norm_multiplier_UI, norm_multiplier_VI, 0, 0, 0,
								 N_local_x, N_local_y, N_local_z, N_incl, N_azimuth, N_frequencies, stride_x, stride_y,
								 stride_z, stride_incl, stride_azimuth, stride_frequencies, stride_stokes);

	write_3d_field_block_to_hdf5(file_id, comm, normalized_output, N_x, N_y, N_z, N_incl, N_azimuth, N_frequencies,
								 N_local_x, N_local_y, N_local_z, local_start_x, local_start_y, local_start_z,
								 &output_field);
}
