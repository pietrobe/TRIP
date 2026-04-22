#include "RT_problem.hpp"
#include "thdf.h"

#include <cstdio>

namespace
{
	//////////////////////////////////////////////////////////////////////////////
	// write_shared_hdf5_metadata
	//////////////////////////////////////////////////////////////////////////////
	void
	write_shared_hdf5_metadata(const hid_t file_id, RT_problem &rt_problem, const int N_x, const int N_y, const int N_z,
							   const int N_incl, const int N_azimuth, const int N_frequencies)
	{
		std::vector<double>		frequencies(rt_problem.nu_grid_.begin(), rt_problem.nu_grid_.end());
		THDF_frequencies_grid_t frequencies_grid;
		frequencies_grid.N_frequencies = N_frequencies;
		frequencies_grid.frequencies   = frequencies.data();

		THDF_write_frequencies_grid_to_hdf5(file_id, &frequencies_grid);

		THDF_angular_grid_t angular_grid;
		angular_grid.N_azimuthal_angles	  = N_azimuth;
		angular_grid.N_inclination_angles = N_incl;
		angular_grid.N_directions		  = N_incl * N_azimuth;

		std::vector<double> theta_grid_rad(rt_problem.theta_grid_.begin(), rt_problem.theta_grid_.end());
		std::vector<double> chi_grid_rad(rt_problem.chi_grid_.begin(), rt_problem.chi_grid_.end());

		std::vector<int> inclinations_indices(angular_grid.N_directions);
		std::vector<int> azimuthal_indices(angular_grid.N_directions);

		for (int i = 0; i < N_incl; i++)
		{
			for (int j = 0; j < N_azimuth; j++)
			{
				inclinations_indices[i * N_azimuth + j] = i;
				azimuthal_indices[i * N_azimuth + j]	= j;
			}
		}

		angular_grid.inclination_angles	  = theta_grid_rad.data();
		angular_grid.azimuthal_angles	  = chi_grid_rad.data();
		angular_grid.inclinations_indices = inclinations_indices.data();
		angular_grid.azimuthal_indices	  = azimuthal_indices.data();
		THDF_write_angular_grid_to_hdf5(file_id, &angular_grid);

		THDF_geometry_3D_t geometry_3D;

		std::vector<double> heights(rt_problem.depth_grid_.begin(), rt_problem.depth_grid_.end());

		geometry_3D.N_x		= N_x;
		geometry_3D.N_y		= N_y;
		geometry_3D.N_z		= N_z;
		geometry_3D.heights = heights.data();
		geometry_3D.delta	= rt_problem.L_;

		THDF_write_geometry_3D_to_hdf5(file_id, &geometry_3D);
	}

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

	/////////////////////////////////////////////////////////////////////////////
	// write_3d_field_block_mpi_conditional
	/////////////////////////////////////////////////////////////////////////////
	template <typename IN_REAL>
	void
	write_3d_field_block_mpi_conditional(hid_t			 file_id,			 //
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
										 IN_REAL		*stokes_IQUI,		 //
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
	{
		if constexpr (sizeof(IN_REAL) == sizeof(THDF_float32_t))
		{
			write_3d_field_block_mpi_f32(
				file_id, comm, normalized_output, N_x, N_y, N_z, N_incl, N_azimuth, N_frequencies, N_stokes, N_local_x,
				N_local_y, N_local_z, local_start_x, local_start_y, local_start_z, stokes_IQUI, stokes_I, stokes_QI,
				stokes_UI, stokes_VI, norm_multiplier_I, norm_multiplier_QI, norm_multiplier_UI, norm_multiplier_VI,
				stride_x, stride_y, stride_z, stride_incl, stride_azimuth, stride_frequencies, stride_stokes);
		}
		else
		{
			write_3d_field_block_mpi(
				file_id, comm, normalized_output, N_x, N_y, N_z, N_incl, N_azimuth, N_frequencies, N_stokes, N_local_x,
				N_local_y, N_local_z, local_start_x, local_start_y, local_start_z, stokes_IQUI, stokes_I, stokes_QI,
				stokes_UI, stokes_VI, norm_multiplier_I, norm_multiplier_QI, norm_multiplier_UI, norm_multiplier_VI,
				stride_x, stride_y, stride_z, stride_incl, stride_azimuth, stride_frequencies, stride_stokes);
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	// copy_3D_block_field
	//////////////////////////////////////////////////////////////////////////////
	template <typename IN_REAL>
	void
	copy_3D_block_field(THDF_3D_field_t &field,					//
						const IN_REAL	*stokes_IQUV,			//
						THDF_float32_t	*stokes_out_I,			//
						THDF_float32_t	*stokes_out_QI,			//
						THDF_float32_t	*stokes_out_UI,			//
						THDF_float32_t	*stokes_out_VI,			//
						const int		 N0_in_local_x,			//
						const int		 N0_in_local_y,			//
						const int		 N0_in_local_z,			//
						const int		 N_local_x,				//
						const int		 N_local_y,				//
						const int		 N_local_z,				//
						const int		 N_incl,				//
						const int		 N_azimuth,				//
						const int		 N_frequencies,			//
						const int		 stride_in_x,			//
						const int		 stride_in_y,			//
						const int		 stride_in_z,			//
						const int		 stride_in_incl,		//
						const int		 stride_in_azimuth,		//
						const int		 stride_in_frequencies, //
						const int		 stride_in_stokes)
	{
		field.norm_multiplier_I	 = nullptr;
		field.norm_multiplier_QI = nullptr;
		field.norm_multiplier_UI = nullptr;
		field.norm_multiplier_VI = nullptr;

		field.stokes_I	= stokes_out_I;
		field.stokes_QI = stokes_out_QI;
		field.stokes_UI = stokes_out_UI;
		field.stokes_VI = stokes_out_VI;

		auto get_idx = [&](int i, int j, int k, int incl, int az, int f, int s)
		{
			return get_input_index(i, j, k, incl, az, f, s, N0_in_local_x, N0_in_local_y, N0_in_local_z, stride_in_x,
								   stride_in_y, stride_in_z, stride_in_incl, stride_in_azimuth, stride_in_frequencies,
								   stride_in_stokes);
		};

		auto get_output_index = [&](const int i,	//
									const int j,	//
									const int k,	//
									const int incl, //
									const int az,	//
									const int f)	//
		{
			return ((((i * N_local_y + j) //
						  * N_local_z
					  + k)
						 * N_incl
					 + incl)
						* N_azimuth
					+ az)
					   * N_frequencies
				   + f;
		};

		for (int i = 0; i < N_local_x; ++i)
		{
			for (int j = 0; j < N_local_y; ++j)
			{
				for (int k = 0; k < N_local_z; ++k)
				{
					for (int incl = 0; incl < N_incl; ++incl)
					{
						for (int az = 0; az < N_azimuth; ++az)
						{
							for (int f = 0; f < N_frequencies; ++f)
							{
								const int out_idx = get_output_index(i, j, k, incl, az, f);

								const int in_I	= get_idx(i, j, k, incl, az, f, 0);
								const int in_QI = get_idx(i, j, k, incl, az, f, 1);
								const int in_UI = get_idx(i, j, k, incl, az, f, 2);
								const int in_VI = get_idx(i, j, k, incl, az, f, 3);

								const THDF_float_t stokes_I_val = stokes_IQUV[in_I];
								const THDF_float_t inv_pc_stokes_I_value =
									(stokes_I_val != 0.0) ? 100.0 / stokes_I_val : 0.0;

								field.stokes_I[out_idx] = static_cast<THDF_float32_t>(stokes_I_val);

								field.stokes_QI[out_idx] = static_cast<THDF_float32_t>(stokes_IQUV[in_QI]		 //
																					   * inv_pc_stokes_I_value); //

								field.stokes_UI[out_idx] = static_cast<THDF_float32_t>(stokes_IQUV[in_UI]		 //
																					   * inv_pc_stokes_I_value); //

								field.stokes_VI[out_idx] = static_cast<THDF_float32_t>(stokes_IQUV[in_VI]		 //
																					   * inv_pc_stokes_I_value); //
							}
						}
					}
				}
			}
		}
	}

} // namespace

////////////////////////////////////////////////////////////////////////
// Main function to write the whole 3D field to HDF5 using the FAL-C input format
////////////////////////////////////////////////////////////////////////
int
write_3D_whole_field_falp_hdf5(RT_problem &rt_problem, const std::string &output_file, //
							   const int step_z_, bool normalized_output)
{
	const double t_start = MPI_Wtime();

	const int step_z = step_z_;

	const int mpi_rank = rt_problem.mpi_rank_;

	{
		// void *_p = malloc(1);
		// if (_p) free(_p);
		// fprintf(stderr, "[HEAP PROBE r%d] write_3D ENTRY\n", mpi_rank);
		// fflush(stderr);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0)
	{
		fprintf(stderr, "[HEAP PROBE r0] write_3D ENTRY barrier passed\n");
		fflush(stderr);
	}

	const auto [i_start, j_start, k_start] = rt_problem.space_grid_->getGhostMargins();

	const int size_i = rt_problem.space_grid_->getLocalSizeX();
	const int size_j = rt_problem.space_grid_->getLocalSizeY();
	const int size_k = rt_problem.space_grid_->getLocalSizeZ();

	const int local_start_x = rt_problem.space_grid_->local_to_global_coordinate(0, i_start);
	const int local_start_y = rt_problem.space_grid_->local_to_global_coordinate(1, j_start);
	const int local_start_z = rt_problem.space_grid_->local_to_global_coordinate(2, k_start);

	const int N_x = rt_problem.N_x_;
	const int N_y = rt_problem.N_y_;
	const int N_z = rt_problem.N_z_;

	const int N_stokes		= 4;
	const int N_frequencies = rt_problem.N_nu_;
	const int N_azimuth		= rt_problem.N_chi_;
	const int N_incl		= rt_problem.N_theta_;

	// const int stride_stokes		 = 1;
	// const int stride_frequencies = N_stokes;
	// const int stride_azimuth	 = N_stokes * N_frequencies;
	// const int stride_incl		 = N_stokes * N_frequencies * N_azimuth;
	// const int stride_z			 = N_stokes * N_frequencies * N_azimuth * N_incl;
	// const int stride_y			 = N_stokes * N_frequencies * N_azimuth * N_incl * size_k;
	// const int stride_x			 = N_stokes * N_frequencies * N_azimuth * N_incl * size_k * size_j;

	int max_size_k;
	int min_size_k;

	// Diagnostic: all ranks report their layout before any allocation.
	// This confirms the function is entered and prints stride info that would reveal a mismatch.
	{
		const int block_sz	= N_stokes * N_incl * N_azimuth * N_frequencies;
		const int stride_z_ = block_sz;
		const int stride_y_ = block_sz * size_k;
		const int stride_x_ = block_sz * size_k * size_j;
		// fprintf(stderr,
		// 		"[DIAG r%d] write_3D ENTER: domain=(%d,%d,%d) gstart=(%d,%d,%d) "
		// 		"N=(stokes=%d freq=%d az=%d incl=%d) block_sz=%d "
		// 		"stride_z=%d stride_y=%d stride_x=%d\n",
		// 		mpi_rank, size_i, size_j, size_k, local_start_x, local_start_y, local_start_z, N_stokes, N_frequencies,
		// 		N_azimuth, N_incl, block_sz, stride_z_, stride_y_, stride_x_);
		fflush(stderr);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0)
	{
		fprintf(stderr, "[DIAG r0] write_3D: all ranks passed ENTER barrier\n");
		fflush(stderr);
	}

	const double total_size = (double)(N_x * N_y * N_z) * (double)(N_incl * N_azimuth * N_frequencies * N_stokes);

	MPI_Allreduce(&size_k, &max_size_k, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&size_k, &min_size_k, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

	if (mpi_rank == 0)
	{
		fprintf(stderr, "[DIAG r0] write_3D: max_size_k=%d min_size_k=%d\n", max_size_k, min_size_k);
		fflush(stderr);
	}

	{ // Collective file creation:
		double t0 = MPI_Wtime();

		if (mpi_rank == 0)
		{
			fprintf(stderr, "[DIAG r0] write_3D: calling THDF_create_zslab_file_single\n");
			fflush(stderr);
		}

		const int flag_create_file = THDF_create_zslab_file_single(MPI_COMM_WORLD, output_file.c_str(), normalized_output,
																   N_x, N_y, N_z, N_incl, N_azimuth, N_frequencies);

		if (flag_create_file < 0)
		{
			fprintf(stderr, "Rank %d: create failed\n", mpi_rank);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank == 0)
		{
			printf("[t] create: %.3f s\n", MPI_Wtime() - t0);
			fflush(stdout);
			fprintf(stderr, "[DIAG r0] write_3D: THDF_create_zslab_file_single OK\n");
			fflush(stderr);
		}
	}

	if (mpi_rank == 0)
	{
		fprintf(stderr, "[DIAG r0] write_3D: serial metadata write start\n");
		fflush(stderr);

		hid_t file_id = H5Fopen(output_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		if (file_id < 0)
		{
			fprintf(stderr, "Rank 0: serial open failed\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		write_shared_hdf5_metadata(file_id, rt_problem,				  //
								   N_x, N_y, N_z,					  //
								   N_incl, N_azimuth, N_frequencies); //

		H5Fflush(file_id, H5F_SCOPE_GLOBAL);
		H5Fclose(file_id);

		fprintf(stderr, "[DIAG r0] write_3D: serial metadata write done\n");
		fflush(stderr);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0)
	{
		fprintf(stderr, "[DIAG r0] write_3D: entering parallel open+write section\n");
		fflush(stderr);
	}

	{
		const double t0		 = MPI_Wtime();
		hid_t		 fapl	 = build_fapl_from_env(MPI_COMM_WORLD);
		hid_t		 file_id = H5Fopen(output_file.c_str(), H5F_ACC_RDWR, fapl);
		H5Pclose(fapl); // fapl is copied into file_id by H5Fopen, safe to close
		if (file_id < 0)
		{
			fprintf(stderr, "Rank %d: H5Fopen failed\n", mpi_rank);
			fflush(stderr);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		if (mpi_rank == 0)
		{
			fprintf(stderr, "[DIAG r0] write_3D: H5Fopen OK, opening zslab handler\n");
			fflush(stderr);
		}

		THDF_zslab_handler_t *handler = THDF_open_zslab_handler_single(file_id,							  //
																	   normalized_output,				  //
																	   N_x, N_y, N_z,					  //
																	   N_incl, N_azimuth, N_frequencies); //
		if (!handler)
		{
			fprintf(stderr, "Rank %d: open_handler failed\n", mpi_rank);
			fflush(stderr);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank == 0)
		{
			printf("[write 3D] open handler: %.3f s\n", MPI_Wtime() - t0);
			fflush(stdout);
			fprintf(stderr, "[DIAG r0] write_3D: zslab handler open OK\n");
			fflush(stderr);
		}

		// Prepare data buffer for writing slices
		std::vector<THDF_float32_t> stokes_I_buffer(size_i * size_j * step_z * N_incl * N_azimuth * N_frequencies);
		std::vector<THDF_float32_t> stokes_Q_buffer(size_i * size_j * step_z * N_incl * N_azimuth * N_frequencies);
		std::vector<THDF_float32_t> stokes_U_buffer(size_i * size_j * step_z * N_incl * N_azimuth * N_frequencies);
		std::vector<THDF_float32_t> stokes_V_buffer(size_i * size_j * step_z * N_incl * N_azimuth * N_frequencies);

		const double wt0 = MPI_Wtime();

		const int N_local_z = size_k;
		const int N_local_x = size_i;
		const int N_local_y = size_j;

		for (int local_zi = 0; local_zi < max_size_k; local_zi += step_z)
		{
			const int is_writer	 = (local_zi < N_local_z);
			const int current_nz = is_writer ? ((local_zi + step_z <= N_local_z) ? step_z : (N_local_z - local_zi)) : 0;
			const int global_start_z = local_start_z + local_zi;

			THDF_3D_field_t field;
			if (is_writer)
			{
				const int stride_stokes		 = 1;
				const int stride_frequencies = N_stokes;
				const int stride_azimuth	 = N_stokes * N_frequencies;
				const int stride_incl		 = N_stokes * N_frequencies * N_azimuth;
				const int stride_z			 = N_stokes * N_frequencies * N_azimuth * N_incl;
				// stride_y/x must span the FULL local z/y extent (size_k, size_j), not just
				// the current slab height — Field::block() lays out as [x][y][z_full][block].
				const int stride_y = N_stokes * N_frequencies * N_azimuth * N_incl * size_k;
				const int stride_x = N_stokes * N_frequencies * N_azimuth * N_incl * size_k * size_j;

				fprintf(stderr,
						"[DIAG r%d] copy slab local_zi=%d current_nz=%d "
						"strides: z=%d y=%d x=%d buf_cap=%d\n",
						mpi_rank, local_zi, current_nz, stride_z, stride_y, stride_x,
						size_i * size_j * step_z * N_incl * N_azimuth * N_frequencies);
				fflush(stderr);

				Real *stokes_IQUV = rt_problem.I_field_->block(0, 0, local_zi);
				copy_3D_block_field<Real>(field,				  //
										  stokes_IQUV,			  // Input data in float 64 bits
										  stokes_I_buffer.data(), //
										  stokes_Q_buffer.data(), //
										  stokes_U_buffer.data(), //
										  stokes_V_buffer.data(), //
										  0,					  //
										  0,					  //
										  0,					  //
										  size_i,				  //
										  size_j,				  //
										  current_nz,			  //
										  N_incl,				  //
										  N_azimuth,			  //
										  N_frequencies,		  //
										  stride_x,				  //
										  stride_y,				  //
										  stride_z,				  //
										  stride_incl,			  //
										  stride_azimuth,		  //
										  stride_frequencies,	  //
										  stride_stokes);		  //

				fprintf(stderr, "[DIAG r%d] copy slab local_zi=%d DONE\n", mpi_rank, local_zi);
				fflush(stderr);
			}

			double tw0 = MPI_Wtime();

			fprintf(stderr, "[DIAG r%d] THDF_write_zslab_block_single local_zi=%d current_nz=%d\n", mpi_rank, local_zi,
					current_nz);
			fflush(stderr);

			/* ALL ranks call this — non-writers pass current_nz=0 */
			THDF_write_zslab_block_single(handler,										//
										  is_writer ? &field : nullptr,					//
										  local_start_x, local_start_y, global_start_z, /* GLOBAL z offset */
										  N_local_x, N_local_y, current_nz, N_incl, N_azimuth, N_frequencies);

			fprintf(stderr, "[DIAG r%d] THDF_write_zslab_block_single local_zi=%d DONE\n", mpi_rank, local_zi);
			fflush(stderr);

			MPI_Barrier(MPI_COMM_WORLD);
			if (mpi_rank == 0)
			{
				printf("[write 3D] z=%d nz=%d: %.3f s\n", global_start_z, current_nz, MPI_Wtime() - tw0);
				const double gb_written = (double)N_x * N_y * current_nz * N_incl * N_azimuth * N_frequencies * N_stokes
										  * sizeof(THDF_float32_t) / (1024.0 * 1024.0 * 1024.0);
				printf("[write 3D] local written: %.3f GB  approx. peak-bandwidth: %.3f GB/s\n", gb_written,
					   gb_written / (MPI_Wtime() - tw0));
				fflush(stdout);
			}
		} // End of loop over z-slices

		MPI_Barrier(MPI_COMM_WORLD);

		if (mpi_rank == 0)
		{
			printf("[write 3D] write loop: %.3f s\n", MPI_Wtime() - wt0);
			double gb_written = (double)total_size * sizeof(THDF_float32_t) / (1024.0 * 1024.0 * 1024.0);
			printf("[write 3D] total written: %.1f GB  effective loop-bandwidth: %.1f GB/s\n", gb_written,
				   gb_written / (MPI_Wtime() - wt0));
			fflush(stdout);
			fprintf(stderr, "[DIAG r0] write_3D: write loop DONE\n");
			fflush(stderr);
		}

		{
			double t0 = MPI_Wtime();
			if (mpi_rank == 0)
			{
				fprintf(stderr, "[DIAG r0] THDF_close_zslab_handler_single\n");
				fflush(stderr);
			}
			THDF_close_zslab_handler_single(handler);
			if (mpi_rank == 0)
			{
				fprintf(stderr, "[DIAG r0] H5Fflush\n");
				fflush(stderr);
			}
			H5Fflush(file_id, H5F_SCOPE_GLOBAL);
			MPI_Barrier(MPI_COMM_WORLD);
			if (mpi_rank == 0)
			{
				printf("[write 3D] H5Fflush: %.3f s\n", MPI_Wtime() - t0);
				fflush(stdout);
			}
			t0 = MPI_Wtime();
			if (mpi_rank == 0)
			{
				fprintf(stderr, "[DIAG r0] H5Fclose\n");
				fflush(stderr);
			}
			H5Fclose(file_id);
			MPI_Barrier(MPI_COMM_WORLD);
			if (mpi_rank == 0)
			{
				printf("[write 3D] H5Fclose: %.3f s\n", MPI_Wtime() - t0);
				fflush(stdout);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	const double t_total = MPI_Wtime() - t_start;
	if (mpi_rank == 0)
	{
		const double gb = total_size * sizeof(THDF_float32_t) / (1024.0 * 1024.0 * 1024.0);
		const int	 Px = rt_problem.mpi_size_x_;
		const int	 Py = rt_problem.mpi_size_y_;
		const int	 Pz = rt_problem.mpi_size_z_;

		printf("[write 3D] total: %.3f s\n", t_total);
		printf("[write 3D] N_x=%d  N_y=%d  N_z=%d  N_incl=%d  N_azimuth=%d  N_frequencies=%d  N_stokes=%d\n", N_x, N_y,
			   N_z, N_incl, N_azimuth, N_frequencies, N_stokes);
		printf("[write 3D] total size: %.1f GB\n", gb);
		printf("[write 3D] max local z-slab: %d  min local z-slab: %d\n", max_size_k, min_size_k);
		printf("[write 3D] decomposition: (%d, %d, %d) for domain (%d, %d, %d)\n", Px, Py, Pz, N_x, N_y, N_z);
		printf("[write 3D] output: %.1f GB  TTS-bandwidth: %.1f GB/s\n", gb, gb / t_total);
	}

	return 0;
}

/////////////////////////////////////////////////////////////////////////
// write_3D_whole_field_hdf5
/////////////////////////////////////////////////////////////////////////
int
write_3D_whole_field_hdf5(RT_problem &rt_problem, const std::string &output_file, //
						  const int step_z_, bool normalized_output)
{
	int step_z = step_z_; // TODO ...

	const int mpi_rank = rt_problem.mpi_rank_;

	const auto [i_start, j_start, k_start] = rt_problem.space_grid_->getGhostMargins();

	const int size_i = rt_problem.space_grid_->getLocalSizeX();
	const int size_j = rt_problem.space_grid_->getLocalSizeY();
	const int size_k = rt_problem.space_grid_->getLocalSizeZ();

	const int local_start_x = rt_problem.space_grid_->local_to_global_coordinate(0, i_start);
	const int local_start_y = rt_problem.space_grid_->local_to_global_coordinate(1, j_start);
	const int local_start_z = rt_problem.space_grid_->local_to_global_coordinate(2, k_start);

	const int N_x = rt_problem.N_x_;
	const int N_y = rt_problem.N_y_;
	const int N_z = rt_problem.N_z_;

	const int N_stokes		= 4;
	const int N_frequencies = rt_problem.N_nu_;
	const int N_azimuth		= rt_problem.N_chi_;
	const int N_incl		= rt_problem.N_theta_;

	const int stride_stokes		 = 1;
	const int stride_frequencies = N_stokes;
	const int stride_azimuth	 = N_stokes * N_frequencies;
	const int stride_incl		 = N_stokes * N_frequencies * N_azimuth;
	const int stride_z			 = N_stokes * N_frequencies * N_azimuth * N_incl;
	const int stride_y			 = N_stokes * N_frequencies * N_azimuth * N_incl * size_k;
	const int stride_x			 = N_stokes * N_frequencies * N_azimuth * N_incl * size_k * size_j;

	int max_size_k;
	int min_size_k;

	MPI_Allreduce(&size_k, &max_size_k, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&size_k, &min_size_k, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

	if (mpi_rank == 0)
	{
		const double hdf_open_and_handler_start = MPI_Wtime();
		// hid_t					 file_id					= THDF_open_file(output_file.c_str());
		if (THDF_create_3D_field_handler_mpi(output_file.c_str(), //
											 normalized_output,	  //
											 N_x, N_y, N_z,		  //
											 N_incl, N_azimuth, N_frequencies)
			< 0)
		{
			std::cerr << "[rank 0] THDF_create_3D_field_handler_mpi failed\n";
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		hid_t file_id = H5Fopen(output_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

		write_shared_hdf5_metadata(file_id,							  //
								   rt_problem,						  //
								   N_x, N_y, N_z,					  //
								   N_incl, N_azimuth, N_frequencies); //

		H5Fclose(file_id);
	} // End of writing shared parameters from MPI rank 0,

	MPI_Barrier(MPI_COMM_WORLD);
	// Write field data in parallel.

	std::vector<THDF_float32_t> slice_data_I_buffer(size_i * size_j * step_z * N_incl * N_azimuth * N_frequencies);
	std::vector<THDF_float32_t> slice_data_Q_buffer(size_i * size_j * step_z * N_incl * N_azimuth * N_frequencies);
	std::vector<THDF_float32_t> slice_data_U_buffer(size_i * size_j * step_z * N_incl * N_azimuth * N_frequencies);
	std::vector<THDF_float32_t> slice_data_V_buffer(size_i * size_j * step_z * N_incl * N_azimuth * N_frequencies);

	if (max_size_k < step_z)
	{
		if (mpi_rank == 0)
		{
			std::cout << "Warning: max local Z size (" << max_size_k << ") is smaller than step_z (" << step_z
					  << "). Adjusting step_z to " << max_size_k << "." << std::endl;
		}
		step_z = max_size_k;
	}

	const double hdf_write_start = MPI_Wtime();

	for (int local_k = 0; local_k < max_size_k; local_k += step_z)
	{
		const bool is_writer	= local_k < size_k ? true : false;
		const int  size_slice_k = local_k + step_z <= size_k ? step_z : size_k - local_k;
		const int  id_z			= local_start_z + local_k;

		MPI_Comm write_communicator;
		MPI_Comm_split(MPI_COMM_WORLD, is_writer ? 1 : MPI_UNDEFINED, mpi_rank, &write_communicator);

		if (is_writer)
		{
			hid_t w_fapl = H5Pcreate(H5P_FILE_ACCESS);
			H5Pset_fapl_mpio(w_fapl, write_communicator, MPI_INFO_NULL);
			hid_t write_file_id = H5Fopen(output_file.c_str(), H5F_ACC_RDWR, w_fapl);
			H5Pclose(w_fapl);

			Real *stokes_IQUI = rt_problem.I_field_->block(0, 0, local_k);

			write_3d_field_block_mpi_conditional(write_file_id,		 //
												 write_communicator, //
												 normalized_output,	 // If true the output will be normalized.
												 N_x,				 // Global sizes of the field in x y z directions
												 N_y,				 // and in inclination, azimuth, frequencies
												 N_z,				 //
												 N_incl,			 //
												 N_azimuth,			 //
												 N_frequencies,		 //
												 N_stokes,			 // ... and number of Stokes parameters (4)
												 size_i,			 // Size of the slice
												 size_j,			 // in x y z directions
												 size_slice_k,		 // size of the current slice in z direction
												 local_start_x,		 // Global start coordinates of the slice
												 local_start_y,		 // in x y z directions
												 id_z,		  // Global start coordinate of the slice in z direction
												 stokes_IQUI, // Pointer to the slice data in memory
												 slice_data_I_buffer.data(), //
												 slice_data_Q_buffer.data(), //
												 slice_data_U_buffer.data(), //
												 slice_data_V_buffer.data(), //
												 nullptr,					 // Pointers to normalization multipliers
												 nullptr,					 // (set to nullptr if not needed)
												 nullptr,					 //
												 nullptr,					 //
												 stride_x,					 // Strides in memory for each dimension
												 stride_y,			 // (x, y, z, incl, azimuth, frequencies, stokes)
												 stride_z,			 //
												 stride_incl,		 //
												 stride_azimuth,	 //
												 stride_frequencies, //
												 stride_stokes);	 //

			H5Fclose(write_file_id);
			MPI_Barrier(write_communicator); // Ensure all writers have finished before we free the communicator
			MPI_Comm_free(&write_communicator);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD); // Ensure all processes have finished writing before we proceed
	const double hdf_write_end = MPI_Wtime();
	if (mpi_rank == 0)
		std::cout << "[rank 0] HDF5 3D write time: " << hdf_write_end - hdf_write_start << " s" << std::endl;

	return 0;
}