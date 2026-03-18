#include "RT_problem.hpp"
#include "thdf.h"

#include <cstdio>

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
		const double			 hdf_open_and_handler_start = MPI_Wtime();
		hid_t					 file_id					= THDF_open_file(output_file.c_str());
		THDF_3D_field_handler_t *field_handler = THDF_create_3D_field_handler_mpi(file_id,							 //
																				  normalized_output,				 //
																				  N_x, N_y, N_z,					 //
																				  N_incl, N_azimuth, N_frequencies); //
		const double			 hdf_open_and_handler_end = MPI_Wtime();
		std::cout << "[rank 0] THDF_write 3D field + THDF_create_3D_field_handler_mpi: "
				  << hdf_open_and_handler_end - hdf_open_and_handler_start << " s" << std::endl;

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

		THDF_close_3D_field_handler_mpi(field_handler);
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
			Real *stokes_IQUI = rt_problem.I_field_->block(0, 0, local_k);

			write_3d_field_block_mpi(output_file.c_str(),		 //
									 write_communicator,		 //
									 normalized_output,			 // If true the output will be normalized.
									 N_x,						 // Global sizes of the field in x y z directions
									 N_y,						 // and in inclination, azimuth, frequencies
									 N_z,						 //
									 N_incl,					 //
									 N_azimuth,					 //
									 N_frequencies,				 //
									 N_stokes,					 // ... and number of Stokes parameters (4)
									 size_i,					 // Size of the slice
									 size_j,					 // in x y z directions
									 size_slice_k,				 // size of the current slice in z direction
									 local_start_x,				 // Global start coordinates of the slice
									 local_start_y,				 // in x y z directions
									 id_z,						 // Global start coordinate of the slice in z direction
									 stokes_IQUI,				 // Pointer to the slice data in memory
									 slice_data_I_buffer.data(), //
									 slice_data_Q_buffer.data(), //
									 slice_data_U_buffer.data(), //
									 slice_data_V_buffer.data(), //
									 nullptr,					 // Pointers to normalization multipliers
									 nullptr,					 // (set to nullptr if not needed)
									 nullptr,					 //
									 nullptr,					 //
									 stride_x,					 // Strides in memory for each dimension
									 stride_y,					 // (x, y, z, incl, azimuth, frequencies, stokes)
									 stride_z,					 //
									 stride_incl,				 //
									 stride_azimuth,			 //
									 stride_frequencies,		 //
									 stride_stokes);			 //

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