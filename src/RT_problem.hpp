#ifndef RT_problem_hpp
#define RT_problem_hpp

#include "GridManager/GridManager.hpp"
#include "RT_utility.hpp"
#include "domain_decomposition_3D.hpp"
#include "hdf_atmos_cpp.hpp"
#include "thdf.h"

// Type for storing fields
// using Real_t = Real;

using Grid_ptr_t  = std::shared_ptr<Grid3D>;
using Field_ptr_t = std::shared_ptr<Field>;

typedef const std::string input_string;

#define GAUSS_TO_LARMOR_FREQUENCY(BB__) ((BB__) * (1399600.0)) // [Gauss] -> [Hz]

inline std::string
emissivity_model_to_string(const emissivity_model_t &model)
{
	switch (model)
	{
		case emissivity_model_t::NONE:
			return "NONE";
			break;
		case emissivity_model_t::CRD_limit_VHP:
		case emissivity_model_t::CRD_limit:
			return "CRD";
			break;
		case emissivity_model_t::PRD:
		case emissivity_model_t::PRD_NORMAL:
		case emissivity_model_t::PRD_FAST:
			return "PRD";
			break;
		case emissivity_model_t::PRD_AA:
		case emissivity_model_t::PRD_AA_MAPV:
			return "PRD_AA";
			break;
		case emissivity_model_t::ZERO:
			return "CONTINUUM";
			break;
		default:
			return "UNKNOWN";
			break;
	}
}

inline std::string
emissivity_model_to_string_long(const emissivity_model_t &model)
{
	switch (model)
	{
		case emissivity_model_t::NONE:
			return "NONE";
			break;
		case emissivity_model_t::CRD_limit:
			return "CRD_limit";
			break;
		case emissivity_model_t::CRD_limit_VHP:
			return "CRD_limit_VHP";
			break;
		case emissivity_model_t::PRD:
			return "PRD";
			break;
		case emissivity_model_t::PRD_NORMAL:
			return "PRD_NORMAL";
			break;
		case emissivity_model_t::PRD_FAST:
			return "PRD_FAST";
			break;
		case emissivity_model_t::PRD_AA:
			return "PRD_AA";
			break;
		case emissivity_model_t::PRD_AA_MAPV:
			return "PRD_AA_MAPV";
			break;
		case emissivity_model_t::ZERO:
			return "CONTINUUM";
			break;
		default:
			return "UNKNOWN";
			break;
	}
}

class RT_problem
{
   public:
	// Deafult constructor
	RT_problem(const AppConfig &cfg)
	{
		Real start = MPI_Wtime();

		// assign MPI varaibles
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

		// set input file
		cfg_ = cfg;

		if (mpi_rank_ == 0 and verbose_)
		{
			std::cout << "\n~~~~~~ MPI size = " << mpi_size_ << " ~~~~~~" << std::endl << std::endl;
			std::cout << "Emissivity Model long:  " << emissivity_model_to_string_long(cfg_.emissivity_model)
					  << std::endl;
			std::cout << "Emissivity Model short: " << emissivity_model_to_string(cfg_.emissivity_model) << std::endl;
		}

		// set flags
		use_PORTA_input_	= not(cfg_.input_directory.string().find("FAL-C") != std::string::npos);
		use_magnetic_field_ = cfg_.use_B;
		emissivity_model_	= cfg_.emissivity_model;
		enable_continuum_	= cfg_.enable_continuum;
		use_1_5D_approx_	= cfg_.use_1_5D_approx;
		use_bulk_velocity_	= cfg_.use_Vb;
		verbose_			= cfg_.verbose;

		// atom quantities
		mass_ = cfg.atom.mass;
		Aul_  = cfg.atom.Aul;
		S2_   = cfg.atom.S2;
		Ll2_  = cfg.atom.Ll2; 
		Lu2_  = cfg.atom.Lu2; 
		
		Eu_vec_ = cfg.atom.Eu_vec;
		El_vec_ = cfg.atom.El_vec;
		gu_vec_ = cfg.atom.gu_vec;
		gl_vec_ = cfg.atom.gl_vec;
		
		use_uniform_magnetic_field_ = cfg_.set_uniform_B;
		if (use_uniform_magnetic_field_)
		{
			uniform_magnetic_field_value_ = cfg_.B_field[0];
			uniform_magnetic_field_theta_ = cfg_.B_field[1];
			uniform_magnetic_field_chi_	  = cfg_.B_field[2];
		}

		// set input files
		const auto input_file_path		    = cfg_.input_directory / cfg_.input_file;
		const auto frequencies_input_path = cfg_.input_directory / cfg_.frequency_file;

		// read input frequencies from separate file
		read_frequency(frequencies_input_path.string().c_str());

		// read atmoshperic data depending on input format
		if (use_PORTA_input_)
		{
			if (cfg_.input_cul.string().empty() || cfg_.input_qel.empty() || cfg_.input_llp.empty())
			{
				if (mpi_rank_ == 0) std::cout << "Using PORTA PMD input file ONLY:  " << input_file_path << std::endl;

				if (hdf_atmos_cpp::is_valid_hdf5_file(input_file_path))
				{
					this->read_3D_h5(input_file_path.string(), this->cfg_, true);
				}
				else
				{
					read_3D(input_file_path.string().c_str());
				}
			}
			else
			{
				if (mpi_rank_ == 0) std::cout << "Using PORTA PMD + CUL + QEL + LLP + BACK input files" << std::endl;

				auto input_cul_path	 = cfg_.input_directory / cfg_.input_cul;
				auto input_qel_path	 = cfg_.input_directory / cfg_.input_qel;
				auto input_llp_path	 = cfg_.input_directory / cfg_.input_llp;
				auto input_back_path = cfg_.input_directory / cfg_.input_back;

				read_3D(input_file_path.string().c_str(), input_cul_path.string().c_str(),
						input_qel_path.string().c_str(), input_llp_path.string().c_str(),
						input_back_path.string().c_str());
			}
		}
		else // FAL-C input
		{
			if (mpi_rank_ == 0) std::cout << "Using FAL-C input format" << std::endl << std::endl;

			read_FAL_C();
		}

		// timing
		if (mpi_rank_ == 0 and verbose_) printf("Reading input time:\t\t%g (seconds)\n", MPI_Wtime() - start);
		start = MPI_Wtime();

		// set up atomic quantities
		set_up_atom();

		// precompute
		set_up();

		if (mpi_rank_ == 0 and verbose_) printf("Set up time:\t\t%g (seconds)\n", MPI_Wtime() - start);

		// print info
		if (verbose_) print_info();
	}

	// DEPRECATED CONSTRUCTORS

	// constructor for PORTA input file with additional inputs
	RT_problem(const char *filename_pmd, const char *filename_cul, const char *filename_qel, const char *filename_llp,
			   const char *filename_back, input_string input_path_frequency,
			   const emissivity_model_t emissivity_model_arg, const bool use_magnetic_field = false)
	{
		Real start = MPI_Wtime();

		// assign MPI varaibles
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

		if (mpi_rank_ == 0)
		{
			std::cout << "\n~~~~~~ MPI size = " << mpi_size_ << " ~~~~~~" << std::endl << std::endl;
			std::cout << "Emissivity Model long:  " << emissivity_model_to_string_long(emissivity_model_arg) << std::endl;
			std::cout << "Emissivity Model short: " << emissivity_model_to_string(emissivity_model_arg) << std::endl;
		}

		// set flags
		use_PORTA_input_	= true;
		use_magnetic_field_ = use_magnetic_field;
		emissivity_model_	= emissivity_model_arg;

		// frequency grid is not contained in PORTA input (but can be computed from T_ref)
		// const bool use_wavelength = false; // TEST
		read_frequency(input_path_frequency);
		read_3D(filename_pmd, filename_cul, filename_qel, filename_llp, filename_back);

		// timing
		if (mpi_rank_ == 0) printf("Reading input time:\t\t%g (seconds)\n", MPI_Wtime() - start);
		start = MPI_Wtime();

		// precompute
		set_up();

		print_info();

		if (mpi_rank_ == 0) printf("Set up time:\t\t%g (seconds)\n", MPI_Wtime() - start);
	}

	/////////////////////////////////////////////////////////////////////////
	// constructor for PORTA input file
	RT_problem(const char *PORTA_input, input_string input_path_frequency, const emissivity_model_t emissivity_model_arg,
			   const bool use_magnetic_field = false)
	{
		Real start = MPI_Wtime();

		// assign MPI varaibles
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

		if (mpi_rank_ == 0) std::cout << "\n~~~~~~ MPI size = " << mpi_size_ << " ~~~~~~" << std::endl;

		// set flags
		use_PORTA_input_	= true;
		use_magnetic_field_ = use_magnetic_field;
		emissivity_model_	= emissivity_model_arg;

		// New method to set the emissivity model
		// emissivity_model emissivity_model_ = emissivity_model::NONE;

		// frequency grid is not contained in PORTA input (but can be computed from T_ref)
		// const bool use_wavelength = false; // TEST
		read_frequency(input_path_frequency);
		read_3D(PORTA_input);

		// timing
		// MPI_Barrier(space_grid_->raw_comm());
		Real end = MPI_Wtime();
		if (mpi_rank_ == 0) printf("Reading input time:\t\t%g (seconds)\n", end - start);
		start = MPI_Wtime();

		// precompute
		set_up();

		print_info();

		// MPI_Barrier(space_grid_->raw_comm());
		end = MPI_Wtime();
		if (mpi_rank_ == 0) printf("Set up time:\t\t%g (seconds)\n", end - start);
	}

	// constructor
	RT_problem(input_string input_path, const int N_theta, const int N_chi,
			   const emissivity_model_t emissivity_model_arg = emissivity_model_t::NONE,
			   const bool				use_magnetic_field	 = false)
	{
		Real start = MPI_Wtime();

		// set emissivity model
		emissivity_model_ = emissivity_model_arg;

		// set CRD flag
		use_magnetic_field_ = use_magnetic_field;

		// assign MPI varaibles
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

		if (mpi_rank_ == 0) std::cout << "\n~~~~~~ MPI size = " << mpi_size_ << " ~~~~~~" << std::endl;
		if (mpi_rank_ == 0) std::cout << "\n=========== Reading input files ===========\n" << std::endl;

		// TODO: now hardcoded (put everything in input_path file)
		// N_x_ = std::sqrt(mpi_size_)/2;
		N_x_ = 10;
		N_y_ = N_x_;
		L_	 = 400.0;
		// const Real L_tot = 1000.0;
		// L_ = L_tot/N_x_;

		// reading some input
		read_atom(input_path + "/atom.dat");
		read_depth(input_path + "/atmosphere.dat");
		read_frequency(input_path + "/frequency.dat");

		// set sizes and directions grids and weigths
		set_theta_chi_grids(N_theta, N_chi);
		set_sizes();

		// init grid
		space_grid_ = std::make_shared<Grid3D>(MPI_COMM_WORLD, N_x_, N_y_, N_z_);
		space_grid_->print_info();

		// init fields
		allocate_fields();

		// init atmospheric quantities
		allocate_atmosphere();

		// read atm data (needs grid object)
		read_atmosphere_1D(input_path + "/atmosphere.dat"); // NOTE: solar surface for space index k = 0
		read_bulk_velocity_1D(input_path + "/bulk_velocity.dat");
		read_magnetic_field_1D(input_path + "/magnetic_field.dat");

		read_continumm_1D(input_path + "/continuum/continuum_scat_opac.dat",
						  input_path + "/continuum/continuum_tot_opac.dat",
						  input_path + "/continuum/continuum_therm_emiss.dat");
		// precompute
		set_up();

		MPI_Barrier(space_grid_->getComm());
		Real end = MPI_Wtime();
		if (mpi_rank_ == 0) printf("Set up time:\t\t%g (seconds)\n", end - start);
	}

	inline void
	read_FAL_C()
	{
		// Set grid horizontal dimensions
		N_x_ = cfg_.N_x;
		N_y_ = cfg_.N_y;
		L_	 = cfg_.L;

		input_string input_path = (cfg_.input_directory / cfg_.input_file).string();

		// reading some input
		read_atom(input_path + "/atom.dat");
		read_depth(input_path + "/atmosphere.dat");

		// set sizes and directions grids and weigths
		set_theta_chi_grids(cfg_.N_theta, cfg_.N_chi);
		set_sizes();

		// create space grid
		set_grid_partition();
		space_grid_ = std::make_shared<Grid3D>(MPI_COMM_WORLD, N_x_, N_y_, N_z_,
											   std::array<PetscInt, 3>{mpi_size_x_, mpi_size_y_, mpi_size_z_});

		print_info();

		// init fields
		allocate_fields();

		// init atmospheric quantities
		allocate_atmosphere();

		// read atm data (needs grid object)
		read_atmosphere_1D(input_path + "/atmosphere.dat"); // NOTE: solar surface for space index k = 0
		read_bulk_velocity_1D(input_path + "/bulk_velocity.dat");
		read_magnetic_field_1D(input_path + "/magnetic_field.dat");

		read_continumm_1D(input_path + "/continuum/continuum_scat_opac.dat",
						  input_path + "/continuum/continuum_tot_opac.dat",
						  input_path + "/continuum/continuum_therm_emiss.dat");
	}

	// convert block index to to local ones = [j_theta, k_chi, n_nu, i_stokes]
	// TODO this will be deleted
	inline std::vector<int>
	block_to_local(const int block_index)
	{
		std::vector<int> local_indeces;
		local_indeces.reserve(4);

		const int i_stokes = block_index % 4;
		const int n		   = (block_index / 4) % N_nu_;
		const int k		   = (block_index / (4 * N_nu_)) % N_chi_;
		const int j		   = (block_index / (4 * N_nu_ * N_chi_)) % N_theta_;

		local_indeces.push_back(j);
		local_indeces.push_back(k);
		local_indeces.push_back(n);
		local_indeces.push_back(i_stokes);

		return local_indeces;
	}

	// convert block index to to local ones = [j_theta, k_chi, n_nu]
	// TODO this will be deleted
	inline std::vector<int>
	block_to_local_unpol(const int block_index)
	{
		std::vector<int> local_indeces;
		local_indeces.reserve(3);

		const int n = block_index % N_nu_;
		const int k = (block_index / N_nu_) % N_chi_;
		const int j = (block_index / (N_nu_ * N_chi_)) % N_theta_;

		local_indeces.push_back(j);
		local_indeces.push_back(k);
		local_indeces.push_back(n);

		return local_indeces;
	}

	// extract horizontal plane k from field
	std::vector<Real>
	extract_plane_k(const Field_ptr_t field, const int k_global);

	// convert local indeces to block one (of fields) for the first Stokes parameter and vice versa
	// TODO this will be deleted
	inline int
	local_to_block(const int j, const int k, const int n)
	{
		return 4 * (N_nu_ * (N_chi_ * j + k) + n);
	}

	// convert local indeces to block one (of fields)
	// TODO this will be deleted
	inline int
	local_to_block_unpol(const int j, const int k, const int n)
	{
		return N_nu_ * (N_chi_ * j + k) + n;
	}


	// set Q,U,V = 0
	inline void set_zero_polarization(Field_ptr_t field)
	{		
    	space_grid_->parallel_for([&](int i, int j, int k) 
    	{         	                        
	        for (int b = 0; b < block_size_; b = b + 4) 
	        {
	            field->block(i,j,k)[b + 1] = 0;                
	            field->block(i,j,k)[b + 2] = 0;                
	            field->block(i,j,k)[b + 3] = 0;                
	        }                                                	        
	    }); 
	}


	// print I_field on surface
	void
	print_surface_profile(const Field_ptr_t field, const int i_stoke = 0, const int i_space = 0, const int j_space = 0,
						  const int j_theta = 0, const int k_chi = 0) const;

	void
	print_surface_QI_profile(const Field_ptr_t field, const int i_space = 0, const int j_space = 0, const int j_theta = 0,
							 const int k_chi = 0, const int i_stokes = 1, const bool center_line = false) const;

	void
	print_surface_QI_point(const int i_space = 0, const int j_space = 0, const int j_theta = 0, const int k_chi = 0,
						   const int n_nu = 0, const int i_stokes = 1) const;

	void
	print_profile(const Field_ptr_t field, const int i_stoke = 0, const int i_space = 0, const int j_space = 0,
				  const int k_space = 0, const int j_theta = 0, const int k_chi = 0) const;

	// write surface profile in file in MATLAB format in onw point and all surface
	void
	write_surface_point_profiles(input_string file_name, const int i_space, const int j_space) const;
	void
	write_surface_point_profiles_Omega(input_string file_name, const int i_space, const int j_space) const;
	void
	write_surface_profiles(input_string file_name) const;
	void
	write_surface_profiles_Omega(input_string file_name) const;

	void
	write_surface_point_profiles_csv(input_string file_name, const int i_space, const int j_space,
									 const unsigned int precision = 14) const;

	void
	write_angular_grid_csv(input_string file_name, const int i_space, const int j_space,
						   const unsigned int precision = 14) const;

	void
	write_frequencies_grid_csv(input_string file_name, const int i_space, const int j_space,
							   const unsigned int precision = 14) const;

	void
	write_surface_point_profiles_Omega_csv(input_string file_name, const int i_space, const int j_space,
										   const unsigned int precision) const;

	std::tuple<int, int>
	make_write_surface_MPI_Comm(const MPI_Comm MPI_Comm_MAIN, MPI_Comm &write_comm);

	int
	accumulate_surface_domain_data(std::vector<Real> &surface_data_I,  //
								   std::vector<Real> &surface_data_Q,  //
								   std::vector<Real> &surface_data_U,  //
								   std::vector<Real> &surface_data_V); //

	/**
	 * Accumulate J, K, Q values over a given subdomain
	 * @param x_strat Starting index in the local x-direction
	 * @param y_strat Starting index in the local y-direction
	 * @param z_strat Starting index in the local z-direction
	 * @param x_end Ending index in the local x-direction (exclusive)
	 * @param y_end Ending index in the local y-direction (exclusive)
	 * @param z_end Ending index in the local z-direction (exclusive)
	 * @param JKQ_real Vector to accumulate real parts of JKQ values
	 * @param JKQ_imag Vector to accumulate imaginary parts of JKQ values
	 * @param KQ_matrix Matrix (9x2) to put KQ matrix values (the number of KQ values is 9, namely the size divided by 2).
	 * @param KQ_matrix_real Matrix (6x2) to put the used KQ values for the real part after the compression.
	 * @param KQ_matrix_imag Matrix (3x2) to put the used KQ values for the imaginary part after the compression.
	 * @note the number of KQ values after the compression can be obtained from the size of KQ_matrix_real and
	 * KQ_matrix_imag divided by 2 (wich is 6 and 3 respectively).
	 * @return MPI status code
	 */
	int
	accumulate_JKQ_values(const int						   x_strat,		   //
						  const int						   y_strat,		   //
						  const int						   z_strat,		   //
						  const int						   x_end,		   //
						  const int						   y_end,		   //
						  const int						   z_end,		   //
						  std::vector<THDF_JKQ_float_t>	  &JKQ_real,	   //
						  std::vector<THDF_JKQ_float_t>	  &JKQ_imag,	   //
						  std::vector<THDF_JKQ_n_float_t> &JKQ_real_norm,  //
						  std::vector<THDF_JKQ_n_float_t> &JKQ_imag_norm); //

	int
	get_KQ_values(std::vector<int> &KQ_values,					//
				  std::vector<int> &KQ_values_real_compressed,	//
				  std::vector<int> &KQ_values_imag_compressed); //

	int
	write_JKQ_field_hdf5(const std::string &output_file); //

	int																	//
	write_angular_frequency_grids_hdf5(const std::string &output_file); //

	int
	write_emergent_field_hdf5(const std::string &output_file, MPI_Comm write_comm, std::vector<Real> &surface_data_I,
							  std::vector<Real> &surface_data_Q, std::vector<Real> &surface_data_U,
							  std::vector<Real> &surface_data_V);

	int
	accumulate_surface_profiles_Omega_domain_data(std::vector<Real> &surface_data_I,	   //
												  std::vector<Real> &surface_data_Q,	   //
												  std::vector<Real> &surface_data_U,	   //
												  std::vector<Real> &surface_data_V);	   //
	int																					   //
	write_beams_frequency_grids_Omega_hdf5(const std::vector<BeamDirection> &beams,		   //
										   const std::string				&output_file); //

	int
	write_emergent_field_Omega_hdf5(const std::string				 &output_file,	   //
									MPI_Comm						  write_comm,	   //
									const std::vector<BeamDirection> &beams,		   //
									const int						  beam_index,	   //
									std::vector<Real>				 &surface_data_I,  //
									std::vector<Real>				 &surface_data_Q,  //
									std::vector<Real>				 &surface_data_U,  //
									std::vector<Real>				 &surface_data_V); //
	// MPI varables
	int mpi_rank_;
	int mpi_size_;

	// procs in each dimension // TODO these will be useless

	domain_decomposition_3D::DomainInfo mpi_decomposition_;

	int mpi_size_x_;
	int mpi_size_y_;
	int mpi_size_z_;

	// flag to enable continuum
	bool enable_continuum_;

	// Use the 1.5 approximation
	bool use_1_5D_approx_;

	// input file
	AppConfig cfg_;

	// print
	bool verbose_;

	// emissivity model
	emissivity_model_t emissivity_model_ = emissivity_model_t::NONE;

	// spatial grid
	Grid_ptr_t space_grid_;

	// auxiliary grids
	std::vector<double> depth_grid_; // ordering: [surf,...,deep] in km
	std::vector<double> nu_grid_;
	std::vector<double> theta_grid_;
	std::vector<double> mu_grid_;
	std::vector<double> chi_grid_;

	// horixontal spacing
	double L_; // L = dx = dy

	// Legendre and trapezoidal weights
	std::vector<double> w_theta_;
	std::vector<double> w_chi_;

	// grids sizes
	int N_x_;
	int N_y_;
	int N_z_;
	int N_s_; // N_s_ = N_x_ * N_y_ * N_z_

	int N_theta_;
	int N_chi_;
	int N_dirs_; // N_dirs_ = N_theta_ * N_chi_;
	int N_nu_;

	int N_pol_ = 4; // ADDED

	int		 local_size_; // == tot_size_ con mpi_size_ = 1
	PetscInt block_size_; // 4 * N_nu_ * N_theta_ * N_chi_;
	PetscInt tot_size_;	  // N_s_ * block_size;

	// unpolarized sizes (normal ones divided by 4) // UNUSED when introducing the new field class, these values could be
	// avoided
	PetscInt local_size_unpolarized_; // == tot_size_ con mpi_size_ = 1
	PetscInt block_size_unpolarized_; // N_nu_ * N_theta_ * N_chi_;
	PetscInt tot_size_unpolarized_;	  // N_s_ * block_size;

	// unknown quantities
	Field_ptr_t I_field_; // intensity
	Field_ptr_t S_field_; // source function

	// PETSc data structure for intensity
	Vec I_vec_; // TODO useless

	// propagation matrix entries
	Field_ptr_t eta_field_;
	Field_ptr_t rho_field_;

	// reduced fields for unpolarized CRD or AA solutions
	Field_ptr_t I_unpol_field_;
	Field_ptr_t S_unpol_field_;

	// atmospheric quantities
	// Field_ptr_t Nu_;   // upper level populations
	Field_ptr_t Nl_;  // lower level populations
	Field_ptr_t T_;	  // temperature
	Field_ptr_t nH_;  // hydrogen number density
	Field_ptr_t xi_;  // microturbulent velocity (a.k.a. non-thermal microscopic velocity)
	Field_ptr_t Cul_; // rate of inelastic de-exciting collisions
	Field_ptr_t Qel_; // rate of elastic collisions // TODO: for memory, maybe this could be removed, leaving only D2_?
	Field_ptr_t a_;	  // damping constant
	Field_ptr_t W_T_; // Wien function

	// magnetic field [nu_L, theta_B_, chi_B_]
	Field_ptr_t B_;

	// bulk velocities [v_b, theta_b, chi_b]
	Field_ptr_t v_b_;

	// quantities depending on position that can be precomputed
	Field_ptr_t Doppler_width_;
	Field_ptr_t k_L_;	  // frequency-integrated absorption coefficient
	Field_ptr_t epsilon_; // thermalization parameter

	// input quantities depending on position and frequency
	Field_ptr_t u_;		// reduced frequencies
	Field_ptr_t sigma_; // continuum coeffs
	Field_ptr_t k_c_;
	Field_ptr_t eps_c_th_;

	// Access to the atomic model parameters
	inline double atomic_mass() const { return mass_; }
	inline double atomic_El()   const { return El_;   }
	inline double atomic_Eu()   const { return Eu_;   }
	inline int    atomic_Jl2()  const {	return Ll2_;  }
	inline int    atomic_Ju2()  const {	return Lu2_;  }
	inline double atomic_gl()   const { return gl_;   } // TODO use gl_vec_
	inline double atomic_gu()   const { return gu_;   } // TODO use gu_vec_
	inline double atomic_Aul()  const { return Aul_;  }
	
	inline Field_ptr_t get_D2() const { return D2_;      }
	inline bool is_two_level()  const { return S2_ == 0; }
	inline int get_S2()         const { return S2_;      }

	inline std::vector<double>	get_El_vec() const	{return El_vec_;}
	inline std::vector<double>	get_Eu_vec() const	{return Eu_vec_;}
	
	inline mdm::md_matrix<double, 1> get_El0() const	{return El0_;}
	inline mdm::md_matrix<double, 1> get_Eu0() const	{return Eu0_;}
	
	bool
	field_is_zero(const Field_ptr_t field);

	// auxiliary_fields for single direction Omega outside theta and chi grids
	void
	allocate_fields_Omega();

	void
	set_eta_and_rhos_Omega(const double theta, const double chi);

	// allocate reduced data structeres
	void
	allocate_unpolarized_fields();
	void
	allocate_J_KQ_field();

	// polarized_to_unpolarized
	void
	polarized_to_unpolarized_field(const Field_ptr_t field, Field_ptr_t field_unpol);

	Field_ptr_t I_field_Omega_; // intensity
	Field_ptr_t S_field_Omega_; // source function
	Field_ptr_t eta_field_Omega_;
	Field_ptr_t rho_field_Omega_;

	inline void
	free_fields_memory()
	{
		if (mpi_rank_ == 0) std::cout << "Freeing RT_problem fields memory..." << std::endl;

		// I_field_.reset();
		// S_field_.reset();
		rho_field_.reset();
		eta_field_.reset();
	}

   private:
	const bool use_ghost_layers_ = false;

	bool use_PORTA_input_	 = false;
	bool use_magnetic_field_ = false;
	bool use_bulk_velocity_	 = false;

	bool use_uniform_magnetic_field_   = false;
	Real uniform_magnetic_field_value_ = 0.0; // Gauss
	Real uniform_magnetic_field_theta_ = 0.0; // rad
	Real uniform_magnetic_field_chi_   = 0.0; // rad

	// physical constants
	const double c_	= 2.99792458e+10;
	const double k_B_ = 1.38065e-16;
	const double h_	= 6.62607e-27;

	// 2-level atom constants
	double mass_;
	double El_ = 0.0;
	double Eu_;
	double gl_;
	double gu_;
	double Aul_;   // Einstein coefficients for spontaneous emission
	double T_ref_; // Reference temperature
	double a_coef_D2_;
	double b_coef_D2_;

	// int  Jl2_;
	// int  Ju2_;
	int Ll2_;
	int Lu2_;
	int S2_;

	std::vector<int>	Jl2_vec_;
	std::vector<int>	Ju2_vec_;
	std::vector<double> gl_vec_;
	std::vector<double> gu_vec_;

	// energy vectors in different formats
	std::vector<double>		  El_vec_;
	std::vector<double>		  Eu_vec_;
	mdm::md_matrix<double, 1> El0_;
	mdm::md_matrix<double, 1> Eu0_;

	// reference frame
	const double gamma_ = 0.5 * M_PI;

	// quantities needed to read 3D PORTA input
	int node_size_;
	int header_size_;

	// atom constant, to precompute
	double nu_0_;

	// depolarizing rate due to elastic collisions
	Field_ptr_t D1_;
	Field_ptr_t D2_;

	// quantities depending on direction
	std::vector<std::vector<std::complex<double>>> T_KQ_; // polarization tensor 

	// allocate grid fields
	void
	allocate_fields();
	void
	allocate_atmosphere();

	// init fields
	void
	init_field(Field_ptr_t input_field, const Real input_value); // UNUSED

	// set grids and sizes
	void
	set_theta_chi_grids(const int N_theta, const int N_chi, const bool double_GL = true);
	void
	set_sizes();

	// menage grid distribution
	void
	set_grid_partition();
	void
	set_3D_decomposition(const int N_x, const int N_y, const int N_z); // UNUSED
	void
	set_3D_decomposition_BLC(const int mpi_rank, const int mpi_size, const int N_x, const int N_y, const int N_z);

	// read inputs
	void
	read_atom(input_string filename);

	// TODO
	void
	set_up_atom();

	void
	read_depth(input_string filename);

	void
	read_frequency(input_string filename, const bool use_wavelength = false);

	void
	read_atmosphere_1D(input_string filename);

	void
	read_bulk_velocity_1D(input_string filename);

	void
	read_magnetic_field_1D(input_string filename);

	void
	read_continumm_1D(input_string filename_sigma, input_string filename_k_c, input_string filename_eps_c_th);

	// read 3D input from pmd file
	void
	read_3D(const char *filename);

	// read 3D input from pmd file
	void
	read_3D_h5(const std::string filename, const AppConfig &cfg, const bool verbose = false);

	void
	read_3D(const char *filename_pmd, const char *filename_cul, const char *filename_qel, const char *filename_llp,
			const char *filename_back);

	std::vector<double>
	read_single_node(MPI_File fh, const int i, const int j, const int k);

	double
	read_single_node_single_field(MPI_File fh, const int i, const int j, const int k);

	void
	read_single_node_triple_field(MPI_File input_file, const int i, const int j, const int k, double &kappa,
								  double &sigma, double &epsilon);

	// compute polarization tensors (vector of six components)
	std::vector<std::complex<double>>
	compute_T_KQ(const int stokes_i, const double theta, const double chi);

	std::complex<double>
	get_TKQi(const int i_stokes, const int K, const int Q, const int j, const int k);

	std::complex<double>
	get_TKQi(const std::vector<std::complex<double>> T_KQ_i, const int K, const int Q);

	void
	set_TKQ_tensor();

	// compute elements of the propagation matrix K
	void
	set_eta_and_rhos();

	// precompute quantities
	void
	set_up();

	// print infos on screen
	void
	print_info() const;
};

#endif
