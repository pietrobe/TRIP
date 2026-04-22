#include "RT_problem.hpp"
#include "hdf_atmos_cpp.hpp"

// read 3D input from pmd file
void
RT_problem::read_3D_h5(const std::string filename, const bool verbose)
{
	if (mpi_rank_ == 0) std::cout << "Reading HDF5 3D data from " << filename << std::endl;

	hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0)
	{
		throw std::runtime_error("Error opening HDF5 file: " + filename);
	}

	//////////////// read atom data ////////////////////////////////

	THDF_atom_two_levels_t atom;

	if (hdf_atmos_cpp::read_atom(file_id, atom) < 0)
	{
		H5Fclose(file_id);
		throw std::runtime_error("Error reading atom data");
	}

	this->gl_	 = atom.g_lower;
	this->gu_	 = atom.g_upper;
	this->El_	 = atom.E_lower;
	this->Eu_	 = atom.E_upper;
	this->Jl2_	 = atom.jl2;
	this->Ju2_	 = atom.ju2;
	this->mass_	 = atom.atomic_mass;
	this->Aul_	 = atom.Aul;
	this->T_ref_ = 5000.0; // reference temperature (never used).

	// Print Atom data for verification
	if (mpi_rank_ == 0 and verbose)
	{
		std::cout << "Atom data read from HDF5 file:" << std::endl;
		std::cout << "  atomic_mass = " << mass_ << std::endl;
		std::cout << "  E_lower =     " << El_ << std::endl;
		std::cout << "  E_upper =     " << Eu_ << std::endl;
		std::cout << "  g_lower =     " << gl_ << std::endl;
		std::cout << "  g_upper =     " << gu_ << std::endl;
		std::cout << "  jl2 =         " << Jl2_ << std::endl;
		std::cout << "  ju2 =         " << Ju2_ << std::endl;
		std::cout << "  Aul =         " << Aul_ << std::endl;
	}

	//////////////// read frequency grid ////////////////////////////////
	// Read geometry and set the sizes:

	if (hdf_atmos_cpp::read_geometry(file_id, this->N_x_, this->N_y_, this->N_z_, this->L_, this->depth_grid_) < 0)
	{
		H5Fclose(file_id);
		throw std::runtime_error("Error reading geometry data");
	}

	this->L_ *= 1e-5; // Convert from cm to km
	std::reverse(this->depth_grid_.begin(), this->depth_grid_.end());

	for (size_t i = 0; i < depth_grid_.size(); ++i) depth_grid_[i] *= 1e-5; // Convert from cm to km

	if (mpi_rank_ == 0 and verbose)
	{
		std::cout << "Geometry read from HDF5 file:" << std::endl;
		std::cout << "  N_x = " << N_x_ << std::endl;
		std::cout << "  N_y = " << N_y_ << std::endl;
		std::cout << "  N_z = " << N_z_ << std::endl;
		std::cout << "  L   = " << L_ << " km" << std::endl;
	}

	/////////////////// set sizes
	this->set_theta_chi_grids(cfg_.N_theta, cfg_.N_chi);
	this->set_sizes();

	// create space grid
	this->set_grid_partition();

	this->space_grid_ =
		std::make_shared<Grid3D>(MPI_COMM_WORLD,																	  //
								 this->N_x_, this->N_y_, this->N_z_,												  //
								 std::array<PetscInt, 3>{{this->mpi_size_x_, this->mpi_size_y_, this->mpi_size_z_}}); //

	// init fields
	this->allocate_fields();

	// init atmospheric quantities
	this->allocate_atmosphere();

	// const int i_start = this->space_grid_->getGlobalStartX();
	// const int j_start = this->space_grid_->getGlobalStartY();
	// const int k_start = this->space_grid_->getGlobalStartZ();

	// const int size_i = this->space_grid_->getLocalSizeX();
	// const int size_j = this->space_grid_->getLocalSizeY();
	// const int size_k = this->space_grid_->getLocalSizeZ();

	// const int k_reverse_start = N_z_ - k_start - size_k;

	// const int read_success = hdf_atmos_cpp::read_atmos_3D_partial_data(file_id, atmos_data,				  //
	// 																   i_start, j_start, k_reverse_start, //
	// 																   size_i, size_j, size_k);			  //

	// const int stride_i = size_j * size_k;
	// const int stride_j = size_k;
	// const int stride_k = 1;

	this->space_grid_ //
		->parallel_for(
			[&](int i, int j, int k)
			{
				// const int dk		 = size_k - 1 - k; // = size_k-1-(k-k_start), descending
				// const int data_index = i * stride_i + j * stride_j + dk;

				const int i_global = space_grid_->local_to_global_coordinate(0, i);
				const int j_global = space_grid_->local_to_global_coordinate(1, j);
				const int k_global = space_grid_->local_to_global_coordinate(2, k);

				const int k_reverse = (N_z_ - k_global - 1);

				std::vector<THDF_atmos_t> atmos_data;

				const int read_success = hdf_atmos_cpp::read_atmos_3D_partial_data(file_id, atmos_data,			  //
																				   i_global, j_global, k_reverse, //
																				   1, 1, 1);					  //

				if (read_success < 0)
				{
					H5Fclose(file_id);
					throw std::runtime_error("Error reading atmospheric data");
				}

				const THDF_atmos_t &data = atmos_data[0];

				// // get global indeces form local ones
				// // space_grid_->getGlobalStartX() + i - space_grid_->getGhostMarginX();
				// const int i_global = this->space_grid_						 //
				// 						 ->local_to_global_coordinate(0, i); //

				// // space_grid_->getGlobalStartY() + j - space_grid_->getGhostMarginY();
				// const int j_global = this->space_grid_						 //
				// 						 ->local_to_global_coordinate(1, j); //

				// // space_grid_->getGlobalStartZ() + k - space_grid_->getGhostMarginZ();
				// const int k_global = this->space_grid_						 //
				// 						 ->local_to_global_coordinate(2, k); //

				if (data.index_i != i_global	  //
					|| data.index_j != j_global	  //
					|| data.index_k != k_reverse) // reversed z is stored in descending order in z

				{
					std::string error_message =
						"Error: Mismatch between expected global indices and those in the HDF5 data. ";
					error_message += "Expected (i, j, k): (" + std::to_string(i_global) + ", " + std::to_string(j_global)
									 + ", " + std::to_string(k_global) + "). ";
					error_message += "Found (index_i, index_j, index_k): (" + std::to_string(data.index_i) + ", "
									 + std::to_string(data.index_j) + ", " + std::to_string(data.index_k) + ").";

					throw std::runtime_error(error_message);
				}

				this->T_->ref(i, j, k)	 = data.temperature;
				this->Nl_->ref(i, j, k)	 = data.pop_lower_level;
				this->a_->ref(i, j, k)	 = data.damping;
				this->Qel_->ref(i, j, k) = data.Qel;
				this->Cul_->ref(i, j, k) = data.Cul;
				xi_->ref(i, j, k)		 = data.vmicro;

				this->epsilon_->ref(i, j, k) = this->Cul_->ref(i, j, k) /				//
											   (this->Cul_->ref(i, j, k) + this->Aul_); //

				// Calculated BUT it should be in the hdf5 file and read directly
				this->D2_->ref(i, j, k) = 0.5 * this->Qel_->ref(i, j, k);

				if (use_magnetic_field_)
				{
					if (use_uniform_magnetic_field_)
					{
						B_->block(i, j, k)[0] = GAUSS_TO_LARMOR_FREQUENCY(uniform_magnetic_field_value_);
						B_->block(i, j, k)[1] = uniform_magnetic_field_theta_;
						B_->block(i, j, k)[2] = uniform_magnetic_field_chi_;
					}
					else
					{
						// convert to spherical coordinates
						auto B_spherical = convert_cartesian_to_spherical(data.magnetic_field_x,  //
																		  data.magnetic_field_y,  //
																		  data.magnetic_field_z); //

						B_->block(i, j, k)[0] =
							GAUSS_TO_LARMOR_FREQUENCY(B_spherical[0]); // converting to Larmor frequency
						B_->block(i, j, k)[1] = B_spherical[1];
						B_->block(i, j, k)[2] = B_spherical[2];
					}
				}
				else
				{
					B_->block(i, j, k)[0] = 0.0;
					B_->block(i, j, k)[1] = 0.0;
					B_->block(i, j, k)[2] = 0.0;
				}

				if (not(this->use_bulk_velocity_))
				{
					v_b_->block(i, j, k)[0] = 0.0;
					v_b_->block(i, j, k)[1] = 0.0;
					v_b_->block(i, j, k)[2] = 0.0;
				}
				else
				{
					// convert to spherical coordinates
					const auto v_spherical	= convert_cartesian_to_spherical(data.bulk_velocity_x,	//
																			 data.bulk_velocity_y,	//
																			 data.bulk_velocity_z); //
					v_b_->block(i, j, k)[0] = v_spherical[0];
					v_b_->block(i, j, k)[1] = v_spherical[1];
					v_b_->block(i, j, k)[2] = v_spherical[2];
				}

				// continuum
				for (int n = 0; n < N_nu_; ++n)
				{
					// hardcoded to 0.0 as in PORTA
					this->sigma_->block(i, j, k)[n]	   = data.c_scat_opacity_sigma_c;
					this->k_c_->block(i, j, k)[n]	   = data.c_tot_opacity_K_c;
					this->eps_c_th_->block(i, j, k)[n] = data.c_therm_emissivity_epsilon_c;
				}
			});

	// Close the HDF5 file
	H5Fclose(file_id);
}