#ifndef RT_problem_hpp
#define RT_problem_hpp

// #include "Utilities.hpp"
// #include "Rotation_matrix.hpp"
#include "sgrid_Core.hpp"

using Real = double;

using Grid_t  = sgrid::Grid<Real, 3>;
using Field_t = sgrid::Field<Grid_t>;

using Grid_ptr_t  = std::shared_ptr<Grid_t>;
using Field_ptr_t = std::shared_ptr<Field_t>;

class RT_problem
{

public:

	// constructor
	RT_problem(const std::string input_path, const size_t N_theta, const size_t N_chi)			   
	{
		// assign MPI varaibles 
    	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
    	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

    	if (mpi_rank_ == 0) std::cout << "\n~~~~~~ MPI size = " << mpi_size_ << " ~~~~~~" << std::endl;		

    	//reading input
    	// if (mpi_rank_ == 0) std::cout << "\n=========== Reading input files ===========\n" << std::endl;				
    	// TODO	

    	// now hardcoded 

    	N_x_ = 100;     
		N_y_ = 100;     
		N_z_ = 100;  
		N_s_ = N_x_ * N_y_ * N_z_;
		

		N_theta_ = N_theta; 
		N_chi_   = N_chi;   
		N_dirs_  = N_theta_ * N_chi_;   
		N_nu_    = 100;

		block_size_ = 4 * N_nu_ * N_theta_ * N_chi_;
		tot_size_   = N_s_ * block_size_;

		Real start = MPI_Wtime();

		// init grid
		space_grid_ = std::make_shared<Grid_t>();
		space_grid_->init(MPI_COMM_WORLD, {(int)N_x_, (int)N_y_, (int)N_z_}, {0, 0, 0});

		// init fields and atmospheric quantities 
		init_fields();	
		init_atmosphere();	

		MPI_Barrier(space_grid_->raw_comm());
	    Real end = MPI_Wtime();
	    Real user_time = end - start;

	    if (mpi_rank_ == 0) {
	    	printf("Device: \"%s\"\n", typeid(DeviceExecutionSpace).name());
      		printf("Setup and kernel call:\t%g (seconds)\n", user_time);
    	}	

    	start = MPI_Wtime();

    	I_field_->exchange_halos();    	
    	S_field_->exchange_halos();    	

	    MPI_Barrier(space_grid_->raw_comm());

	    end = MPI_Wtime();
	    user_time = end - start;

	    if (mpi_rank_ == 0) {
	      printf("Exchange:\t\t%g (seconds)\n", user_time);
	    }   

	}

private:

	// MPI varables
	int mpi_rank_;
	int mpi_size_;	

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

		// depolarizing rate due to elastic collisions
	std::vector<Real> D1_, D2_; 

	// atom constant, to precompute
	Real nu_0_;	

	// spatial grid
	Grid_ptr_t space_grid_; 

	// auxiliary grids	
	std::vector<Real> nu_grid_; 
	std::vector<Real> theta_grid_; 
	std::vector<Real> mu_grid_; 	
	std::vector<Real> chi_grid_; 

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
	Field_ptr_t Nu_;   // upper level populations 
	Field_ptr_t T_;    // temperature 
	Field_ptr_t xi_;   // microturbulent velocity (a.k.a. non-thermal microscopic velocity)		
	Field_ptr_t nu_L_; // Larmor frequency
	Field_ptr_t Cul_;  // rate of inelastic de-exciting collisions
	Field_ptr_t Qel_;  // rate of elastic collisions 
	Field_ptr_t a_;    // damping constant 
	Field_ptr_t W_T_;  // Wien function

	// magnetic field direction, in polar coordinates   
	Field_ptr_t theta_B_; 
	Field_ptr_t chi_B_;   

	// bulk velocities, in polar coordinates
	Field_ptr_t v_b_;       
	Field_ptr_t theta_v_b_; 
	Field_ptr_t chi_v_b_;   
	
	// quantities depending on position that can be precomputed
	Field_ptr_t Doppler_width_;
	Field_ptr_t k_L_;      // frequency-integrated absorption coefficient
	Field_ptr_t epsilon_;  // thermalization parameter

	// input quantities depending on position and frequency 
	Field_ptr_t u_;      // reduced frequencies 
	Field_ptr_t sigma_;  // continuum coeffs
	Field_ptr_t k_c_;
	Field_ptr_t eps_c_th_;

	// quantities depending on position and direction
	Field_ptr_t T_KQ_; // polarization tensor 

	void init_fields();
	void init_atmosphere();
};

#endif 
