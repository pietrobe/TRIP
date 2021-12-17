#ifndef RT_problem_hpp
#define RT_problem_hpp

#include "Utilities.hpp"
#include "sgrid_Core.hpp"

using Grid_t  = sgrid::Grid<Real, 3>; 
using Field_t = sgrid::Field<Grid_t>;

using Grid_ptr_t  = std::shared_ptr<Grid_t>;
using Field_ptr_t = std::shared_ptr<Field_t>;

typedef const std::string input_string;

class RT_problem
{

public:

	// constructor
	RT_problem(input_string input_path, const size_t N_theta, const size_t N_chi)			   
	{
		// assign MPI varaibles 
    	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
    	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

    	if (mpi_rank_ == 0) std::cout << "\n~~~~~~ MPI size = " << mpi_size_ << " ~~~~~~" << std::endl;		
    	if (mpi_rank_ == 0) std::cout << "\n=========== Reading input files ===========\n" << std::endl;				

    	Real start = MPI_Wtime();

    	// TODO: now hardcoded
    	L_   = 1.0;
    	N_x_ = 4;
    	N_y_ = 4;

    	// reading some input
    	read_depth(    input_path + "/atmosphere.dat");	
    	read_frequency(input_path + "/frequency.dat");	
    	    	
    	// set sizes and directions grids and weigths
		set_theta_chi_grids(N_theta, N_chi);
		set_sizes();
								
		// init grid
		space_grid_ = std::make_shared<Grid_t>();
		space_grid_->init(MPI_COMM_WORLD, {(int)N_x_, (int)N_y_, (int)N_z_}, {1, 1, 0});
		
		MPI_Barrier(space_grid_->raw_comm());
		Real end = MPI_Wtime();
	    Real user_time = end - start;

	    if (mpi_rank_ == 0) printf("Init grid:\t%g (seconds)\n", user_time);

	    start = MPI_Wtime();

		// init fields
		allocate_fields();					

		// init atmospheric quantities 
		allocate_atmosphere();	

		// read atm data (needs grid object)
		read_atmosphere_1D(    input_path + "/atmosphere.dat");			 // NOTE: solar surface for space index k = 0
		read_bulk_velocity_1D( input_path + "/bulk_velocity.dat");	
		read_magnetic_field_1D(input_path + "/magnetic_field.dat");
		
		read_continumm_1D(  input_path + "/continuum/continuum_scat_opac.dat", 
						    input_path + "/continuum/continuum_tot_opac.dat",
						    input_path + "/continuum/continuum_therm_emiss.dat");

		MPI_Barrier(space_grid_->raw_comm());
	    end = MPI_Wtime();
	    user_time = end - start;
	    if (mpi_rank_ == 0) printf("Allocate:\t%g (seconds)\n", user_time);
		start = MPI_Wtime();

		// precompute
		set_up();

		MPI_Barrier(space_grid_->raw_comm());
	    end = MPI_Wtime();
	    user_time = end - start;
	    if (mpi_rank_ == 0) printf("Set up:\t%g (seconds)\n", user_time);
    	start = MPI_Wtime();

    	I_field_->exchange_halos();    	
    	S_field_->exchange_halos();    	

	    MPI_Barrier(space_grid_->raw_comm());

	    end = MPI_Wtime();
	    user_time = end - start;

	    if (mpi_rank_ == 0) printf("Exchange:\t\t%g (seconds)\n", user_time);	      

	    print_info();

		// start = MPI_Wtime();

		// T_->write("T.raw");		
		// k_c_->write("k_c_.raw");				

		// MPI_Barrier(space_grid_->raw_comm());
		// end = MPI_Wtime();
		// user_time = end - start;

		// if (mpi_rank_ == 0) printf("Output Time:\t\t%g (seconds)\n", user_time);		   
	}

	// convert block index to to local ones = [j_theta, k_chi, n_nu, i_stokes]
	inline std::vector<size_t> block_to_local(const size_t block_index)
	{
		std::vector<size_t> local_indeces;

		const size_t i_stokes =  block_index % 4;
		const size_t n        = (block_index / 4) % N_nu_;
		const size_t k        = (block_index / (4 * N_nu_)) % N_chi_;
		const size_t j        = (block_index / (4 * N_nu_ * N_chi_)) % N_theta_;
				
		local_indeces.push_back(j);
		local_indeces.push_back(k);
		local_indeces.push_back(n);
		local_indeces.push_back(i_stokes);

		return local_indeces;	
	}

	// convert local indeces to block one (of fields) for the first Stokes parameter and vice versa
	inline size_t local_to_block(const size_t j, const size_t k, const size_t n) { return 4 * ( N_nu_ * ( N_chi_ * j + k ) + n); }
		
	// MPI varables
	int mpi_rank_;
	int mpi_size_;	

	// flag to enable continuum 
	bool enable_continuum_ = true;
	
	// spatial grid
	Grid_ptr_t space_grid_; 	

	// auxiliary grids	
	std::vector<Real> depth_grid_; // in Km 
	std::vector<Real> nu_grid_; 
	std::vector<Real> theta_grid_; 
	std::vector<Real> mu_grid_; 	
	std::vector<Real> chi_grid_; 
	
	// horixontal spacing
	Real L_; // L = dx = dy

	// Legendre and trapezoidal weights 
	std::vector<Real> w_theta_;
	std::vector<Real> w_chi_;

	// grids sizes	
	size_t N_x_;     
	size_t N_y_;     
	size_t N_z_;  
	size_t N_s_; // N_s_ = N_x_ * N_y_ * N_z_

	size_t N_theta_; 
	size_t N_chi_;   
	size_t N_dirs_; // N_dirs_ = N_theta_ * N_chi_;   
	size_t N_nu_;        
	
	size_t block_size_; // 4 * N_nu_ * N_theta_ * N_chi_;
	size_t tot_size_;   // N_s_ * block_size;

	// unknown quantities 
	Field_ptr_t I_field_; // intensity 
	Field_ptr_t S_field_; // source function
	
	// propagation matrix entries 
	Field_ptr_t eta_field_; 
	Field_ptr_t rho_field_;

	// atmospheric quantities
	Field_ptr_t Nl_;   // lower level populations 
	// Field_ptr_t Nu_;   // upper level populations 
	Field_ptr_t T_;    // temperature 
	Field_ptr_t xi_;   // microturbulent velocity (a.k.a. non-thermal microscopic velocity)			
	Field_ptr_t Cul_;  // rate of inelastic de-exciting collisions
	Field_ptr_t Qel_;  // rate of elastic collisions 
	Field_ptr_t a_;    // damping constant 
	Field_ptr_t W_T_;  // Wien function

	// magnetic field [nu_L, theta_B_, chi_B_]
	Field_ptr_t B_; 	

	// bulk velocities [v_b, theta_b, chi_b]
	Field_ptr_t v_b_;       
	
	// quantities depending on position that can be precomputed
	Field_ptr_t Doppler_width_;
	Field_ptr_t k_L_;      // frequency-integrated absorption coefficient
	Field_ptr_t epsilon_;  // thermalization parameter

	// input quantities depending on position and frequency 
	Field_ptr_t u_;      // reduced frequencies 
	Field_ptr_t sigma_;  // continuum coeffs
	Field_ptr_t k_c_;
	Field_ptr_t eps_c_th_;

private:

	// physical and math constants 
	const Real c_   = 2.99792458e+10;
	const Real k_B_ = 1.38065e-16;
	const Real h_   = 6.62607e-27;
	const Real pi_  = 3.1415926535897932384626;	  

	// 2-level atom constants
	const Real mass_ = 40.078;
	const Real El_   = 0.0;
	const Real Eu_   = 23652.304;
	const int Jl_    = 0.0;
	const int Ju_    = 1.0;
	const int gl_    = 0.0;
	const int gu_    = 1.0;
	const Real Aul_  = 2.18e+8;	// Einstein coefficients for spontaneous emission

	// reference frame
	const Real gamma_ = 0.5 * pi_;	  

	// atom constant, to precompute
	Real nu_0_;	

	// depolarizing rate due to elastic collisions
	Field_ptr_t	D1_; 
	Field_ptr_t	D2_;

	// quantities depending on direction
	std::vector<std::vector<std::complex<Real> > > T_KQ_; // polarization tensor 

	// allocate grid fields 
	void allocate_fields();
	void allocate_atmosphere();

	// init fields 
	void init_field(Field_ptr_t input_field, const Real input_value); // TODO remove?
	
	// set grids and sizes
	void set_theta_chi_grids(const size_t N_theta, const size_t N_chi, const bool double_GL = true);
	void set_sizes();

	// read inputs
	void read_depth(            input_string filename);
	void read_frequency(        input_string filename);
	void read_atmosphere_1D(    input_string filename);
	void read_bulk_velocity_1D( input_string filename);	
	void read_magnetic_field_1D(input_string filename);
	
	void read_continumm_1D(input_string filename_sigma, 
						   input_string filename_k_c, 
						   input_string filename_eps_c_th);

	// compute polarization tensors (vector of six components)
	std::vector<std::complex<Real> > compute_T_KQ(const size_t stokes_i, const Real theta, const Real chi);
	std::complex<Real> get_TKQi(const size_t i_stokes, const int K, const int Q, const size_t j, const size_t k);

	// compute elements of the propagation matrix K
	void set_eta_and_rhos();
	
	// precompute quantities
	void set_up();
	
	// print infos on screen
	void const print_info();
};

#endif 
