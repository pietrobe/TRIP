#include <rii_JKQ.h>
#include "RT_problem.hpp"
#include "thdf.h"

//////////////////////////////////////////////////////////////////////////
// Write emergent field to HDF5
// write_emergent_field_hdf5
//////////////////////////////////////////////////////////////////////////
int
write_emergent_field_hdf5(RT_problem &rt_problem, const std::string &output_file)
{
	MPI_Comm write_comm;
	const auto [color, mpi_status] = rt_problem.make_write_surface_MPI_Comm(MPI_COMM_WORLD, write_comm);

	// Note: color == 1 means that the rank is a surface writer
	//       color == 0 means that the rank is NOT a surface writer
	if (color == 0)
	{
		MPI_Barrier(write_comm);
		MPI_Comm_free(&write_comm);
		return EXIT_SUCCESS;
	}

	rt_problem.write_angular_frequency_grids_hdf5(output_file);
	MPI_Barrier(write_comm);

	const int N_x = rt_problem.N_x_;
	const int N_y = rt_problem.N_y_;

	std::vector<double> surface_data_I;
	std::vector<double> surface_data_Q;
	std::vector<double> surface_data_U;
	std::vector<double> surface_data_V;

	const int N_nu	  = rt_problem.N_nu_;
	const int N_theta = rt_problem.N_theta_;
	const int N_chi	  = rt_problem.N_chi_;

	const auto f_dev = rt_problem.I_field_->view_device();
	const auto g_dev = rt_problem.space_grid_->view_device();

	// indeces
	const int i_start = g_dev.margin[0];
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
	const int k_end = k_start + g_dev.dim[2];

	const int size_i = g_dev.dim[0];
	const int size_j = g_dev.dim[1];

	rt_problem.accumulate_surface_domain_data(surface_data_I, surface_data_Q, surface_data_U, surface_data_V);

	int status = rt_problem.write_emergent_field_hdf5(output_file, write_comm, surface_data_I, surface_data_Q,
													  surface_data_U, surface_data_V);

	MPI_Barrier(write_comm);
	MPI_Comm_free(&write_comm);

	return status;
}

/////////////////////////////////////////////////////////////////////////
// Create MPI communicator for writing surface data
// make_write_surface_MPI_Comm
//////////////////////////////////////////////////////////////////////////
std::tuple<int, int>
RT_problem::make_write_surface_MPI_Comm(const MPI_Comm MPI_Comm_MAIN, MPI_Comm &write_comm)
{
	const auto g_dev = space_grid_->view_device();

	const int color = (g_dev.global_coord(2, g_dev.margin[2]) == 0) ? 1 : 0;

	const int mpi_status =						//
		MPI_Comm_split(MPI_Comm_MAIN,			//
					   color,					//
					   mpi_rank_, &write_comm); //

	return std::make_tuple(color, mpi_status);
}

//////////////////////////////////////////////////////////////////////////
// Write emergent angular and frequency grids to HDF5
// write_angular_frequency_grids_hdf5
///////////////////////////////////////////////////////////////////////////
int																			   //
RT_problem::write_angular_frequency_grids_hdf5(const std::string &output_file) //
{
	const auto g_dev = space_grid_->view_device();

	// Only MPI rank 0 AND the process owning global spatial coords (0,0) should create the file
	if (not(this->mpi_rank_ == 0 and						//
			g_dev.global_coord(0, g_dev.margin[0]) == 0 and //
			g_dev.global_coord(1, g_dev.margin[1]) == 0 and //
			g_dev.global_coord(2, g_dev.margin[2]) == 0))	//
	{
		return EXIT_SUCCESS;
	}

	THDF_angular_grid_t angular_grid;

	angular_grid.N_directions		  = N_theta_ / 2 * N_chi_;
	angular_grid.N_inclination_angles = N_theta_ / 2;
	angular_grid.N_azimuthal_angles	  = N_chi_;

	std::vector<double> theta_grid_rad(N_theta_ / 2);
	std::vector<int>	inclinations_indices;
	std::vector<int>	azimuthal_indices(N_chi_);

	for (int j = N_theta_ / 2; j < N_theta_; ++j)
	{
		theta_grid_rad[j - N_theta_ / 2] = theta_grid_[j];
		inclinations_indices.push_back(j);
	}

	angular_grid.inclination_angles = theta_grid_rad.data();
	angular_grid.azimuthal_angles	= chi_grid_.data();

	std::iota(azimuthal_indices.begin(), azimuthal_indices.end(), 0);
	angular_grid.inclinations_indices = inclinations_indices.data();
	angular_grid.azimuthal_indices	  = azimuthal_indices.data();

	hid_t file_id, plist_id;
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	file_id	 = H5Fcreate(output_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

	if (THDF_write_angular_grid_to_hdf5(file_id, &angular_grid) != 0)
	{
		fprintf(stderr, "Error writing emergent angular grid to HDF5 file, rank %d, file: %s:%d\n", mpi_rank_, __FILE__,
				__LINE__);
		return EXIT_FAILURE;
	}

	THDF_frequencies_grid_t freq_grid;
	freq_grid.N_frequencies = N_nu_;
	freq_grid.frequencies	= nu_grid_.data();

	if (THDF_write_frequencies_grid_to_hdf5(file_id, &freq_grid) != 0)
	{
		fprintf(stderr, "Error writing frequencies grid to HDF5 file, rank %d, file: %s:%d\n", mpi_rank_, __FILE__,
				__LINE__);
		return EXIT_FAILURE;
	}

	H5Fclose(file_id);

	return EXIT_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////
// Write emergent field to HDF5
// write_emergent_field_hdf5
//////////////////////////////////////////////////////////////////////////
int
RT_problem::write_emergent_field_hdf5(const std::string &output_file, MPI_Comm write_comm,
									  std::vector<double> &surface_data_I, std::vector<double> &surface_data_Q,
									  std::vector<double> &surface_data_U, std::vector<double> &surface_data_V)
{
	const auto g_dev   = space_grid_->view_device();
	const int  k_start = g_dev.margin[2];

	if (g_dev.global_coord(2, k_start) != 0) return EXIT_SUCCESS;

	const int i_size = g_dev.dim[0];
	const int j_size = g_dev.dim[1];

	const int i_start = g_dev.margin[0];
	const int j_start = g_dev.margin[1];

	const int i_start_global = g_dev.global_coord(0, i_start);
	const int j_start_global = g_dev.global_coord(1, j_start);

	// Create HDF5 file with MPI I/O access property list
	hid_t file_id, plist_id;
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, write_comm, MPI_INFO_NULL);
	file_id = H5Fopen(output_file.c_str(), H5F_ACC_RDWR, plist_id);
	H5Pclose(plist_id);

	THDF_field_handler_t *output_dset_handler = THDF_create_field_handler_mpi(file_id,						//
																			  N_x_, N_y_,					//
																			  N_theta_ / 2, N_chi_, N_nu_); //

	if (output_dset_handler == NULL)
	{
		fprintf(stderr, "Rank %d: Error creating MPI output field dataset\n", mpi_rank_);
		return EXIT_FAILURE;
	}

	THDF_field_t output_field  = {};
	output_field.index_i	   = i_start;
	output_field.index_j	   = j_start;
	output_field.index_incl	   = 0;
	output_field.index_azimuth = 0;

	// assign data pointers
	output_field.stokes_I  = surface_data_I.data();
	output_field.stokes_QI = surface_data_Q.data();
	output_field.stokes_UI = surface_data_U.data();
	output_field.stokes_VI = surface_data_V.data();

	// Write local output field surface data
	THDF_write_field_dataset_to_hdf5(output_dset_handler,				   //
									 &output_field,						   //
									 i_start_global, j_start_global, 0, 0, //
									 i_size, j_size,					   //
									 N_theta_ / 2, N_chi_, N_nu_);		   //

	THDF_close_field_handler_mpi(output_dset_handler);
	H5Fclose(file_id);

	return EXIT_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////
// Accumulate surface data from local domain
// accumulate_surface_data
//////////////////////////////////////////////////////////////////////////
int
RT_problem::accumulate_surface_domain_data(std::vector<double> &surface_data_I, std::vector<double> &surface_data_Q,
										   std::vector<double> &surface_data_U, std::vector<double> &surface_data_V)
{
	const auto f_dev = I_field_->view_device();
	const auto g_dev = space_grid_->view_device();

	// indeces
	const int i_start = g_dev.margin[0];
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];

	const int size_i = g_dev.dim[0];
	const int size_j = g_dev.dim[1];

	int i_global, j_global;

	if (g_dev.global_coord(2, k_start) != 0) return 0;

	const int N_nu	  = this->N_nu_;
	const int N_theta = this->N_theta_;
	const int N_chi	  = this->N_chi_;

	surface_data_I.reserve(N_nu * (N_theta / 2) * N_chi * size_i * size_j);
	surface_data_Q.reserve(N_nu * (N_theta / 2) * N_chi * size_i * size_j);
	surface_data_U.reserve(N_nu * (N_theta / 2) * N_chi * size_i * size_j);
	surface_data_V.reserve(N_nu * (N_theta / 2) * N_chi * size_i * size_j);

	std::cout << "Rank " << mpi_rank_ << " accumulating surface data from domain: i[" << i_start << ", " << i_end
			  << ")  j[" << j_start << ", " << j_end << ")" << std::endl
			  << " global coords: i[" << g_dev.global_coord(0, i_start) << ", " << g_dev.global_coord(0, i_end - 1)
			  << "]  j[" << g_dev.global_coord(1, j_start) << ", " << g_dev.global_coord(1, j_end - 1) << "]"
			  << std::endl;

	for (int i = i_start; i < i_end; ++i)
	{
		// i_global = g_dev.global_coord(0, i);
		for (int j = j_start; j < j_end; ++j)
		{
			// j_global = g_dev.global_coord(1, j);
			for (int j_theta = N_theta_ / 2; j_theta < N_theta_; ++j_theta)
			{
				for (int k_chi = 0; k_chi < N_chi_; ++k_chi)
				{
					const int b_start = local_to_block(j_theta, k_chi, 0);

					for (int b = 0; b < 4 * N_nu_; b = b + 4)
					{
						const double I = f_dev.block(i, j, k_start)[b_start + b + 0];
						const double Q = f_dev.block(i, j, k_start)[b_start + b + 1];
						const double U = f_dev.block(i, j, k_start)[b_start + b + 2];
						const double V = f_dev.block(i, j, k_start)[b_start + b + 3];

						surface_data_I.push_back(I);
						surface_data_Q.push_back(Q / I * 100.0);
						surface_data_U.push_back(U / I * 100.0);
						surface_data_V.push_back(V / I * 100.0);
					}
				}
			}
		}
	}

	return 0;
}

int
RT_problem::accumulate_surface_profiles_Omega_domain_data(std::vector<double> &surface_data_I, //
														  std::vector<double> &surface_data_Q, //
														  std::vector<double> &surface_data_U, //
														  std::vector<double> &surface_data_V)
{
	const auto f_dev = I_field_->view_device();
	const auto g_dev = space_grid_->view_device();
	const auto I_dev = I_field_Omega_->view_device();

	// indeces
	const int i_start = g_dev.margin[0];
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];

	const int size_i = g_dev.dim[0];
	const int size_j = g_dev.dim[1];

	int i_global, j_global;

	if (g_dev.global_coord(2, k_start) != 0) return 0;

	const int N_nu	  = this->N_nu_;
	const int N_theta = this->N_theta_;
	const int N_chi	  = this->N_chi_;

	surface_data_I.reserve(N_nu * size_i * size_j);
	surface_data_Q.reserve(N_nu * size_i * size_j);
	surface_data_U.reserve(N_nu * size_i * size_j);
	surface_data_V.reserve(N_nu * size_i * size_j);

	for (int i = i_start; i < i_end; ++i)
	{
		i_global = g_dev.global_coord(0, i);
		for (int j = j_start; j < j_end; ++j)
		{
			j_global = g_dev.global_coord(1, j);

			for (int b = 0; b < 4 * N_nu_; b = b + 4)
			{
				I_dev.block(i, j, k_start)[b];

				const double I = I_dev.block(i, j, k_start)[b + 0];
				const double Q = I_dev.block(i, j, k_start)[b + 1];
				const double U = I_dev.block(i, j, k_start)[b + 2];
				const double V = I_dev.block(i, j, k_start)[b + 3];

				surface_data_I.push_back(I);
				surface_data_Q.push_back(Q / I * 100.0);
				surface_data_U.push_back(U / I * 100.0);
				surface_data_V.push_back(V / I * 100.0);
			}
		}
	}

	return 0;
}

int																							//
RT_problem::write_beams_frequency_grids_Omega_hdf5(const std::vector<BeamDirection> &beams, //
												   const std::string				&output_file)
{
	const int N_directions = static_cast<int>(beams.size());

	// Only MPI rank 0 AND the process owning global spatial coords (0,0) should create the file
	const auto g_dev = space_grid_->view_device();
	if (not(this->mpi_rank_ == 0 and						//
			g_dev.global_coord(0, g_dev.margin[0]) == 0 and //
			g_dev.global_coord(1, g_dev.margin[1]) == 0 and //
			g_dev.global_coord(2, g_dev.margin[2]) == 0))	//
	{
		return EXIT_SUCCESS;
	}

	std::vector<double> mu_beams(N_directions);
	std::vector<double> chi_beams(N_directions);

	for (int b = 0; b < N_directions; ++b)
	{
		mu_beams[b]	 = beams[b].mu;
		chi_beams[b] = beams[b].chi;
	}

	hid_t file_id, plist_id;
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	file_id	 = H5Fcreate(output_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

	H5Pclose(plist_id);

	THDF_ard_directions_t ard_directions;
	ard_directions.N_directions = N_directions;
	ard_directions.mu			= mu_beams.data();
	ard_directions.chi			= chi_beams.data();

	if (THDF_write_ard_directions(file_id, &ard_directions) != 0)
	{
		fprintf(stderr, "Error writing beams grid to HDF5 file, rank %d, file: %s:%d\n", mpi_rank_, __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	THDF_ard_frequencies_t ard_freqs;

	ard_freqs.N_frequencies = N_nu_;
	ard_freqs.frequencies	= nu_grid_.data();

	if (THDF_write_ard_frequencies(file_id, &ard_freqs) != 0)
	{
		fprintf(stderr, "Error writing ARD frequencies to HDF5 file, rank %d, file: %s:%d\n", mpi_rank_, __FILE__,
				__LINE__);
		return EXIT_FAILURE;
	}

	H5Fclose(file_id);

	return EXIT_SUCCESS;
}

int
RT_problem::write_emergent_field_Omega_hdf5(const std::string				 &output_file,	  //
											MPI_Comm						  write_comm,	  //
											const std::vector<BeamDirection> &beams,		  //
											const int						  beam_index,	  //
											std::vector<double>				 &surface_data_I, //
											std::vector<double>				 &surface_data_Q, //
											std::vector<double>				 &surface_data_U, //
											std::vector<double>				 &surface_data_V)
{
	const auto g_dev = space_grid_->view_device();

	// indeces
	const int i_start = g_dev.margin[0];
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];

	const int size_i = g_dev.dim[0];
	const int size_j = g_dev.dim[1];

	if (g_dev.global_coord(2, k_start) != 0) return EXIT_SUCCESS;

	const int i_global = g_dev.global_coord(0, i_start);
	const int j_global = g_dev.global_coord(1, j_start);

	hid_t file_id, plist_id;
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, write_comm, MPI_INFO_NULL);
	file_id = H5Fopen(output_file.c_str(), H5F_ACC_RDWR, plist_id);

	H5Pclose(plist_id);

	// Create ARD field handler for this direction
	THDF_ard_field_handler *dset_handler =						 //
		THDF_create_ard_field_handler_mpi(file_id,				 //
										  beams[beam_index].mu,	 //
										  beams[beam_index].chi, //
										  this->N_x_,			 //
										  this->N_y_,			 //
										  beams.size(),			 //
										  this->N_nu_);			 //

	if (dset_handler == NULL)
	{
		fprintf(stderr, "Rank %d: Error creating ARD output field dataset\n", mpi_rank_);
		return EXIT_FAILURE;
	}

	THDF_ard_field_t output_field = {};
	output_field.index_i		  = i_global;
	output_field.index_j		  = j_global;
	output_field.index_k		  = 0;
	output_field.beam_index		  = beam_index;
	output_field.N_frequencies	  = this->N_nu_;
	// assign data pointers
	output_field.stokes_I  = surface_data_I.data();
	output_field.stokes_QI = surface_data_Q.data();
	output_field.stokes_UI = surface_data_U.data();
	output_field.stokes_VI = surface_data_V.data();

	// Write local output field surface data
	const int ret = THDF_write_ard_field_to_hdf5(dset_handler,		 //
												 &output_field,		 //
												 i_global, j_global, //
												 size_i, size_j,	 //
												 this->N_nu_);		 //

	if (ret != 0)
	{
		fprintf(stderr, "Rank %d: Error writing ARD output field dataset\n", mpi_rank_);
		return EXIT_FAILURE;
	}

	THDF_close_ard_field_handler_mpi(dset_handler);
	H5Fclose(file_id);

	return EXIT_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////
// Accumulate JKQ values over a given subdomain
// accumulate_JKQ_values
//////////////////////////////////////////////////////////////////////////
int
RT_problem::accumulate_JKQ_values(const int						 x_strat,  //
								  const int						 y_strat,  //
								  const int						 z_strat,  //
								  const int						 x_end,	   //
								  const int						 y_end,	   //
								  const int						 z_end,	   //
								  std::vector<THDF_JKQ_float_t> &JKQ_real, //
								  std::vector<THDF_JKQ_float_t> &JKQ_imag)
{
	const auto f_dev			 = I_field_->view_device();
	auto	   Doppler_width_dev = Doppler_width_->view_device();

	std::vector<double> u_vec;
	u_vec.resize(this->N_nu_);
	for (int i = x_strat; i < x_end; ++i)
	{
		for (int j = y_strat; j < y_end; ++j)
		{
			for (int k = z_strat; k < z_end; ++k)
			{
				const Real_t *block_ptr = f_dev.block(i, j, k);
				const double  dnd		= Doppler_width_dev.ref(i, j, k);

				for (int ui = 0; ui < this->N_nu_; ++ui)
				{
					const double nu = this->nu_grid_[ui];
					u_vec[ui]		= (this->nu_0_ - nu) / dnd;
				}

				auto JKQ_matrix_sh_ptr =												   //
					rii_include::make_JKQ_matrix_norm_comp(block_ptr,					   //
														   u_vec.data(),				   //
														   this->N_nu_,					   //
														   this->N_theta_,				   //
														   this->N_chi_,				   //
														   4,							   //
														   4 * this->N_chi_ * this->N_nu_, //
														   4 * this->N_nu_,				   //
														   1);							   //

				const int size_JKQ_real_comp = JKQ_matrix_sh_ptr->size_KQ_real_compressed();
				const int size_JKQ_imag_comp = JKQ_matrix_sh_ptr->size_KQ_imag_compressed();

				for (int idx = 0; idx < size_JKQ_real_comp; ++idx)
				{
					for (int ui = 0; ui < this->N_nu_; ++ui)
					{
						JKQ_real.push_back(JKQ_matrix_sh_ptr->real_compressed(idx, ui));
					}
				}

				for (int idx = 0; idx < size_JKQ_imag_comp; ++idx)
				{
					for (int ui = 0; ui < this->N_nu_; ++ui)
					{
						JKQ_imag.push_back(JKQ_matrix_sh_ptr->imag_compressed(idx, ui));
					}
				}
			}
		}
	}

	return 0;
}

int
RT_problem::get_KQ_values(std::vector<int> &KQ_values,				   //
						  std::vector<int> &KQ_values_real_compressed, //
						  std::vector<int> &KQ_values_imag_compressed)
{
	rii_include::make_KQ_values(KQ_values,					//
								KQ_values_real_compressed,	//
								KQ_values_imag_compressed); //

	return 0;
}

int
RT_problem::write_JKQ_field_hdf5(const std::string &output_file)
{
	const auto g_dev = space_grid_->view_device();

	// indeces
	const int i_start = g_dev.margin[0];
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
	const int k_end = k_start + g_dev.dim[2];

	const int size_i = g_dev.dim[0];
	const int size_j = g_dev.dim[1];
	const int size_k = g_dev.dim[2];

	std::vector<int> KQ_values;
	std::vector<int> KQ_values_real_compressed;
	std::vector<int> KQ_values_imag_compressed;

	this->get_KQ_values(KQ_values, KQ_values_real_compressed, KQ_values_imag_compressed);

	if (this->mpi_rank_ == 0)
	{
		std::cout << " Writing JKQ field to HDF5 file: " << output_file << std::endl;

		hid_t file_id, plist_id;
		plist_id = H5Pcreate(H5P_FILE_ACCESS);
		file_id	 = H5Fcreate(output_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
		H5Pclose(plist_id);

		THDF_frequencies_grid_t freq_grid;
		freq_grid.N_frequencies = this->N_nu_;
		freq_grid.frequencies	= this->nu_grid_.data();
		if (THDF_write_frequencies_grid_to_hdf5(file_id, &freq_grid) != 0)
		{
			fprintf(stderr, "Error writing frequencies grid to HDF5 file\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		// write JKQ KQ matrix
		THDF_KQ_matrix_t KQ_matrix;
		KQ_matrix.KQ_size				  = KQ_values.size() / 2;
		KQ_matrix.KQ_matrix				  = KQ_values.data();
		KQ_matrix.KQ_compressed_size_real = KQ_values_real_compressed.size() / 2;
		KQ_matrix.KQ_compressed_size_imag = KQ_values_imag_compressed.size() / 2;
		KQ_matrix.KQ_compressed_real	  = KQ_values_real_compressed.data();
		KQ_matrix.KQ_compressed_imag	  = KQ_values_imag_compressed.data();

		if (THDF_write_KQ_matrix_to_hdf5(file_id, &KQ_matrix) != 0)
		{
			fprintf(stderr, "Error writing KQ matrix to HDF5 file\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		H5Fclose(file_id);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	std::vector<THDF_JKQ_float_t> JKQ_real;
	std::vector<THDF_JKQ_float_t> JKQ_imag;

	JKQ_real.reserve(this->N_nu_ * size_i * size_j * size_k * KQ_values_real_compressed.size() / 2);
	JKQ_imag.reserve(this->N_nu_ * size_i * size_j * size_k * KQ_values_imag_compressed.size() / 2);

	this->accumulate_JKQ_values(i_start, j_start, k_start, //
								i_end, j_end, k_end,	   //
								JKQ_real, JKQ_imag);	   //
}
