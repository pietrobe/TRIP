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

	// constructor for PORTA input file
	RT_problem(const char* PORTA_input, input_string input_path_frequency)
	{				
		Real start = MPI_Wtime();

		// assign MPI varaibles 
    	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
    	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

    	if (mpi_rank_ == 0) std::cout << "\n~~~~~~ MPI size = " << mpi_size_ << " ~~~~~~" << std::endl;		

		// set flags    	
    	use_PORTA_input_ = true;
    	use_CRD_limit_   = true;

    	// frequency grid is not conteined in PORTA input (but can be computed from T_ref)
    	read_frequency(input_path_frequency + "/frequency.dat");
    	read_3D(PORTA_input);

    	// timing
    	MPI_Barrier(space_grid_->raw_comm()); Real end = MPI_Wtime(); 	    
	    if (mpi_rank_ == 0) printf("Reading input time:\t\t%g (seconds)\n", end - start);	      		
    	start = MPI_Wtime();

		// precompute
		set_up();

		print_info();	    
	
	   	MPI_Barrier(space_grid_->raw_comm()); end = MPI_Wtime(); 	    
	    if (mpi_rank_ == 0) printf("Set up time:\t\t%g (seconds)\n", end - start);	      		
	}

	// constructor
	RT_problem(input_string input_path, const size_t N_theta, const size_t N_chi)			   
	{
		Real start = MPI_Wtime();

		// assign MPI varaibles 
    	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
    	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

    	if (mpi_rank_ == 0) std::cout << "\n~~~~~~ MPI size = " << mpi_size_ << " ~~~~~~" << std::endl;		
    	if (mpi_rank_ == 0) std::cout << "\n=========== Reading input files ===========\n" << std::endl;				
    
    	// TODO: now hardcoded (put everything in input_path file)    	
    	// N_x_ = std::sqrt(mpi_size_)/2;
    	N_x_ = 1;
    	N_y_ = N_x_;
    	L_   = 100.0;
    	// const double L_tot = 1000.0;
    	// L_ = L_tot/N_x_;

    	// reading some input
    	read_atom(     input_path + "/atom.dat");
    	read_depth(    input_path + "/atmosphere.dat");	
    	read_frequency(input_path + "/frequency.dat");	
    	    	
    	// set sizes and directions grids and weigths
		set_theta_chi_grids(N_theta, N_chi);
		set_sizes();
								
		// init grid
		space_grid_ = std::make_shared<Grid_t>();

		// menage grid distribution // TODO: now bit hardcoded // necessary?
		set_grid_partition();
		space_grid_->init(MPI_COMM_WORLD, {(int)N_x_, (int)N_y_, (int)N_z_}, {1, 1, 0},
									 {mpi_size_x_, mpi_size_y_, mpi_size_z_}, use_ghost_layers_); 
		
		// space_grid_->init(MPI_COMM_WORLD, {(int)N_x_, (int)N_y_, (int)N_z_}, {1, 1, 0}, {}, use_ghost_layers); 

		// init fields
		allocate_fields();				

		// init atmospheric quantities 
		allocate_atmosphere();	

		// read atm data (needs grid object)
		read_atmosphere_1D(    input_path + "/atmosphere.dat");    // NOTE: solar surface for space index k = 0
		read_bulk_velocity_1D( input_path + "/bulk_velocity.dat");	
		read_magnetic_field_1D(input_path + "/magnetic_field.dat");
		
		read_continumm_1D(input_path + "/continuum/continuum_scat_opac.dat", 
						  input_path + "/continuum/continuum_tot_opac.dat",
						  input_path + "/continuum/continuum_therm_emiss.dat");		
		// precompute
		set_up();
		
	    print_info();	    

	   	MPI_Barrier(space_grid_->raw_comm()); Real end = MPI_Wtime(); 	    
	    if (mpi_rank_ == 0) printf("Set up time:\t\t%g (seconds)\n", end - start);	      		
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

	inline void print_PETSc_mem()
	{
		PetscLogDouble space;
   		PetscMemoryGetCurrentUsage(&space);   

   		const double byte_to_GB = 1.0 / (1000 * 1024 * 1024);

    	if (mpi_rank_ == 0) std::cout << "Memory used by PETSc = " << mpi_size_ * byte_to_GB * space << " GB" <<  std::endl;   
	}

	// print I_field on surface 
	void const print_surface_profile(const Field_ptr_t field, 
									 const int i_stoke = 0, const int i_space = 0, const int j_space = 0, 
									 const int j_theta = 0, const int k_chi = 0);


	void const print_surface_QI_profile(const Field_ptr_t field, 
									    const int i_space = 0, const int j_space = 0, 
									    const int j_theta = 0, const int k_chi = 0, const int i_stokes = 1);


	void const print_surface_QI_point(const int i_space = 0, const int j_space = 0, 
								      const int j_theta = 0, const int k_chi = 0, const int n_nu = 0, const int i_stokes = 1);


	void const print_profile(const Field_ptr_t field, 
							 const int i_stoke = 0, const int i_space = 0, const int j_space = 0, const int k_space = 0,
							 const int j_theta = 0, const int k_chi = 0);
		
	// MPI varables
	int mpi_rank_;
	int mpi_size_;	

	// procs in each dimension
	int mpi_size_x_;
	int mpi_size_y_;
	int mpi_size_z_;
	
	// flag to enable continuum 
	bool enable_continuum_ = true;
	
	// flag to use CRD
	bool use_CRD_limit_ = false;		  

	// spatial grid
	Grid_ptr_t space_grid_; 	

	// auxiliary grids	
	std::vector<Real> depth_grid_; // ordering: [surf,...,deep] in km
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
	size_t local_size_; // == tot_size_ con mpi_size_ = 1

	// unknown quantities 
	Field_ptr_t I_field_; // intensity     
	Field_ptr_t S_field_; // source function

	// PETSc data structure for intensity	
	Vec I_vec_;
		
	// propagation matrix entries 
	Field_ptr_t eta_field_; 
	Field_ptr_t rho_field_;

	// atmospheric quantities
	Field_ptr_t Nl_;   // lower level populations 
	// Field_ptr_t Nu_;   // upper level populations 
	Field_ptr_t T_;    // temperature 
	Field_ptr_t xi_;   // microturbulent velocity (a.k.a. non-thermal microscopic velocity)			
	Field_ptr_t Cul_;  // rate of inelastic de-exciting collisions
	Field_ptr_t Qel_;  // rate of elastic collisions // TODO: for memory, maybe this could be removed, leaving only D2_?
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

	// Access to the atomic model parameters
	inline double atomic_mass() const { return mass_;}
	inline double atomic_El()   const { return El_;  }
	inline double atomic_Eu()   const { return Eu_;  }
	inline int    atomic_Jl2()  const { return Jl2_; }
	inline int    atomic_Ju2()  const { return Ju2_; }
	inline double atomic_gl()   const { return gl_;  }
	inline double atomic_gu()   const { return gu_;  } 
	inline double atomic_Aul()  const { return Aul_; }
	inline Field_ptr_t get_D2() const { return D2_;  }

	bool field_is_zero(const Field_ptr_t field);

private:

	const bool use_ghost_layers_ = false;
		  bool use_PORTA_input_  = false;

	// physical constants 
	const Real c_   = 2.99792458e+10;
	const Real k_B_ = 1.38065e-16;
	const Real h_   = 6.62607e-27;

	// 2-level atom constants
	double mass_;
	double El_ = 0.0;
	double Eu_;
	double gl_; 
	double gu_;
	double Aul_;   // Einstein coefficients for spontaneous emission
	double T_ref_; // Reference temperature 
	int Jl2_;
	int Ju2_;
	
	// reference frame
	const Real gamma_ = 0.5 * PI;	  

	// quantities needed to read 3D PORTA input
	int node_size_;
	int header_size_;

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

	// menage grid distribution
	void set_grid_partition();

	// read inputs
	void read_atom(             input_string filename);
	void read_depth(            input_string filename);
	void read_frequency(        input_string filename);
	void read_atmosphere_1D(    input_string filename);
	void read_bulk_velocity_1D( input_string filename);	
	void read_magnetic_field_1D(input_string filename);
	
	void read_continumm_1D(input_string filename_sigma, 
						   input_string filename_k_c, 
						   input_string filename_eps_c_th);

	// read 3D input from pmd file 
	void read_3D(const char* filename);
	std::vector<Real> read_single_node(FILE *f1, const int i, const int j, const int k);

	// compute polarization tensors (vector of six components)
	std::vector<std::complex<Real> > compute_T_KQ(const size_t stokes_i, const Real theta, const Real chi);
	std::complex<Real> get_TKQi(const size_t i_stokes, const int K, const int Q, const size_t j, const size_t k);

	void set_TKQ_tensor();

	// compute elements of the propagation matrix K
	void set_eta_and_rhos();
	
	// precompute quantities
	void set_up();

	// free memory when not needed anymore ----> TODO? ask Simone
	void clean();

	// print infos on screen
	void const print_info();	
};

#endif 
