#include "RT_problem.hpp"
#include "thdf.h"

int
write_3D_whole_field_hdf5(RT_problem &rt_problem, const std::string &output_file, const int step_z_ = 1);

int
write_3D_whole_field_hdf5(RT_problem &rt_problem, const std::string &output_file, const int step_z_)
{
	const int step_z = 1; // TODO ...

	const int mpi_rank = rt_problem.mpi_rank_;

	const auto [i_start, j_start, k_start] = rt_problem.space_grid_->getGhostMargins();

	const int size_i = rt_problem.space_grid_->getLocalSizeX();
	const int size_j = rt_problem.space_grid_->getLocalSizeY();
	const int size_k = rt_problem.space_grid_->getLocalSizeZ();

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

	MPI_Allreduce(&size_k, &max_size_k, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	for (int local_k = 0; local_k < max_size_k; local_k += step_z)
	{
		const bool is_writer	= local_k < size_k ? true : false;
		const int  size_slice_k = step_z;

		MPI_Comm write_communicator;
		MPI_Comm_split(MPI_COMM_WORLD, is_writer ? 1 : MPI_UNDEFINED, mpi_rank, &write_communicator);

		if (is_writer)
		{
			MPI_Barrier(write_communicator); // Ensure all writers have finished before we free the communicator
			MPI_Comm_free(&write_communicator);
		}
	}

	return 0;
}