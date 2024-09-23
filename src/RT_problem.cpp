#include "RT_problem.hpp"

// some bounds for the pmd file
#define PMD_MAIN_HEADER1    89
#define PMD_MAIN_HEADER2    201743
#define TWOLEVEL_HEADER1    4
#define TWOLEVEL_HEADER2    32 

#define GAUSS_TO_LARMOR_FREQUENCY(BB__) ((BB__) * 1399600.0) // [Gauss] -> [Hz]

// read atom and grid quantities 
void RT_problem::read_3D(const char* filename_pmd, const char* filename_cul, const char* filename_qel, const char* filename_llp, const char* filename_back)
{
	if (mpi_rank_ == 0) std::cout << "Reading PORTA input..."   << std::endl;
	if (mpi_rank_ == 0) std::cout << "Reading .pmd input from " << filename_pmd << std::endl;
	if (mpi_rank_ == 0) std::cout << "Reading .cul input from " << filename_cul << std::endl;
	if (mpi_rank_ == 0) std::cout << "Reading .qel input from " << filename_qel << std::endl;
	if (mpi_rank_ == 0) std::cout << "Reading .llp input from " << filename_llp << std::endl;

	const bool zero_velocities = false;
	if (mpi_rank_ == 0 and zero_velocities) std::cout << "WARNING: ZERO velocities HARDCODED!" << std::endl;

	// reading atom and grids from file
	MPI_File fh, f_cul, f_qel, f_llp, f_back;
	MPI_CHECK(MPI_File_open(MPI_COMM_WORLD, filename_pmd, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh));
	MPI_CHECK(MPI_File_open(MPI_COMM_WORLD, filename_cul, MPI_MODE_RDONLY, MPI_INFO_NULL, &f_cul));
	MPI_CHECK(MPI_File_open(MPI_COMM_WORLD, filename_qel, MPI_MODE_RDONLY, MPI_INFO_NULL, &f_qel));
	MPI_CHECK(MPI_File_open(MPI_COMM_WORLD, filename_llp, MPI_MODE_RDONLY, MPI_INFO_NULL, &f_llp));

	// buffers
	Real entry;
	char c;

    // some irrelevant data
	for (int i = 0; i < PMD_MAIN_HEADER1 - 48; i++) { MPI_CHECK(MPI_File_read_all(fh, &c, 1, MPI_CHAR, MPI_STATUS_IGNORE)); }

	// TODO there are not used 
	double Lx, Ly, Lz, x_origin, y_origin, z_origin; // 48 in previuos line for these
	
	MPI_CHECK(MPI_File_read_all(fh, &Lx, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &Ly, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &Lz, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));

	MPI_CHECK(MPI_File_read_all(fh, &x_origin, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &y_origin, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &z_origin, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));

	MPI_CHECK(MPI_File_read_all(fh, &N_x_, 1, MPI_INT, MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &N_y_, 1, MPI_INT, MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &N_z_, 1, MPI_INT, MPI_STATUS_IGNORE));
	
	// WARNING: valid for Lx = Ly = constant 
	L_ = 1e-5 * Lx/N_x_; // conversion from cm to km

	// some irrelevant data (x,y coordinates)	
	int skip_size = 16384 * sizeof(double);
	MPI_CHECK(MPI_File_seek(fh, skip_size, MPI_SEEK_CUR));

	// get z grid (depth)
	for (int i = 0; i < N_z_; ++i)
	{		
		MPI_CHECK(MPI_File_read_all(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
		 
		entry *= 1e-5; // conversion from cm to km			
		depth_grid_.push_back(entry);      		
	}
	
	// reverse to get correct ordering
	reverse(depth_grid_.begin(), depth_grid_.end());

	// some irrelevant data (extra z coordinates)
	skip_size = (8192 - N_z_) * sizeof(double);	
	MPI_CHECK(MPI_File_seek(fh, skip_size, MPI_SEEK_CUR));

	// inclinations and azimuths per octant 	
	MPI_CHECK(MPI_File_read_all(fh, &N_theta_, 1, MPI_INT , MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &N_chi_  , 1, MPI_INT , MPI_STATUS_IGNORE));

	// convert from octant to total 
	N_theta_ *= 2;
	N_chi_   *= 4;

	if (mpi_rank_ == 0) std::cout << "N_theta = " << N_theta_ << ", N_chi = " << N_chi_ << ", from PORTA input." << std::endl;
		
	// some irrelevant data	
	skip_size = PMD_MAIN_HEADER2 - 196608 - 16;
	MPI_CHECK(MPI_File_seek(fh, skip_size, MPI_SEEK_CUR));

    int module_head_size;

	// Size of each grid node	
	MPI_CHECK(MPI_File_read_all(fh, &module_head_size, 1, MPI_INT, MPI_STATUS_IGNORE));

	// Size of each grid node	
	MPI_CHECK(MPI_File_read_all(fh, &node_size_, 1, MPI_INT, MPI_STATUS_IGNORE));

    // Jump to data
    header_size_ = PMD_MAIN_HEADER1 + PMD_MAIN_HEADER2 + 12 + module_head_size;

	// some irrelevant data    
    MPI_CHECK(MPI_File_seek(fh, TWOLEVEL_HEADER1, MPI_SEEK_CUR));

	// reading atomic data    
    MPI_CHECK(MPI_File_read_all(fh, &mass_,  1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &Aul_,   1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &Eu_,    1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &Jl2_,   1, MPI_INT,    MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &Ju2_,   1, MPI_INT,    MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &gl_,    1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &gu_,    1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &T_ref_, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));

    // // change Eu units to [cm-1]     
    // Eu_ /= c_ * h_; 
    
    // hardcoded from FAL-C 
    Eu_ = 23652.304;    
    
  	// some irrelevant data double temp[ny][nx];    matrix of ground (iz=0) for Planckian boundary
	skip_size = (N_x_ * N_y_) * sizeof(double);
	MPI_CHECK(MPI_File_seek(fh, skip_size, MPI_SEEK_CUR));		

	// set angualr grids and sizes and print
	set_theta_chi_grids(N_theta_, N_chi_);
	set_sizes();

	/////////////////// set sizes
	N_nu_ = nu_grid_.size();
	block_size_ = (PetscInt) 4 * N_nu_ * N_theta_ * N_chi_;
	tot_size_   = (PetscInt) N_s_ * block_size_;	
	
	// create space grid
	space_grid_ = std::make_shared<Grid_t>();
	
	// menage grid distribution and init
	set_grid_partition();	
	space_grid_->init(MPI_COMM_WORLD, {N_x_, N_y_, N_z_}, {1, 1, 0},
									 {mpi_size_x_, mpi_size_y_, mpi_size_z_}, use_ghost_layers_); 

	// // TEST
	// space_grid_->init(MPI_COMM_WORLD, {N_x_, N_y_, N_z_}, {1, 1, 0}, {}, use_ghost_layers_); 
		
	// init fields
	allocate_fields();				

	// init atmospheric quantities 
	allocate_atmosphere();	

	// load devices for parallel data structures
	auto T_dev   = T_   ->view_device();
	auto a_dev   = a_   ->view_device();
	auto Nl_dev  = Nl_  ->view_device();
	auto Cul_dev = Cul_ ->view_device();
	auto B_dev   = B_   ->view_device();
	auto v_b_dev = v_b_ ->view_device();
	auto xi_dev  = xi_  ->view_device();
	auto Qel_dev = Qel_ ->view_device();
	auto D2_dev  = D2_  ->view_device();

	auto sigma_dev    = sigma_    ->view_device();
	auto k_c_dev      = k_c_      ->view_device();
	auto eps_c_th_dev = eps_c_th_ ->view_device();
	auto epsilon_dev  = epsilon_  ->view_device();

	auto Doppler_width_dev = Doppler_width_->view_device();

	auto g_dev = space_grid_->view_device();

	if (mpi_rank_ == 0) std::cout << "WARNING: sigma continuum = 0 HARDCODED!" << std::endl;

	// fill field 
	sgrid::parallel_for("READ-ATM1D", space_grid_->md_range(), SGRID_LAMBDA(int i, int j, int k) {

		// get global indeces form local ones
		const int i_global  = g_dev.start[0] + i - g_dev.margin[0];
		const int j_global  = g_dev.start[1] + j - g_dev.margin[1];
		const int k_global  = g_dev.start[2] + k - g_dev.margin[2];
		
		// reversing z index because of input ordering 
		const int k_reverse = (N_z_ - k_global - 1);  

		auto tmp_vector = read_single_node(fh,i_global,j_global,k_reverse);			

		// epsilon_dev.ref(i,j,k) = tmp_vector[0];		
		// Cul_dev.ref(i,j,k)     = tmp_vector[1];		
		T_dev.ref(i,j,k)       = tmp_vector[2];			
		// Nl_dev.ref(i,j,k)      = tmp_vector[9];		
		// a_dev.ref(i,j,k)       = tmp_vector[10];		
		D2_dev.ref(i,j,k)      = tmp_vector[11];

		// hardcoding xi to zero
		xi_dev.ref(i,j,k) = 0; 
		
		// compute Qel and Cul
		Cul_dev.ref(i,j,k) = read_single_node_single_field(f_cul,i_global,j_global,k_reverse);					
		Qel_dev.ref(i,j,k) = read_single_node_single_field(f_qel,i_global,j_global,k_reverse);		
		Nl_dev.ref(i,j,k)  = read_single_node_single_field(f_llp,i_global,j_global,k_reverse);							
		
		// convert to spherical coordinates
		auto B_spherical = convert_cartesian_to_spherical(tmp_vector[3], 
														  tmp_vector[4],
														  tmp_vector[5]);	
		if (use_magnetic_field_)
		{								
			B_dev.block(i, j, k)[0] = B_spherical[0] * 1399600.0; // converting to Larmor frequency					
			B_dev.block(i, j, k)[1] = B_spherical[1]; 					
			B_dev.block(i, j, k)[2] = B_spherical[2]; 

			// // /*  hardcoded B field */ ////////////////////
			/*
			const double B_field_hardcoded = 1000.0; // [Gauss]
			const double theta_B_field = 1.5707963268; // [rad]
			const double chi_B_field = 0.0; // [rad]
			if ( mpi_rank_ == 0 and i == 0 and j == 0 and k == 0) std::cout << "WARNING: HARDCODED B FIELD: of:    " <<  B_field_hardcoded << " Gauss" << std::endl;
			if ( mpi_rank_ == 0 and i == 0 and j == 0 and k == 0) std::cout << "WARNING: HARDCODED B FIELD, theta: " <<  theta_B_field << " rad" << std::endl;
			if ( mpi_rank_ == 0 and i == 0 and j == 0 and k == 0) std::cout << "WARNING: HARDCODED B FIELD, chi:   " <<  chi_B_field << " rad" << std::endl;

			B_dev.block(i, j, k)[0] = GAUSS_TO_LARMOR_FREQUENCY(B_field_hardcoded) ; // converting to Larmor frequency					
			B_dev.block(i, j, k)[1] = theta_B_field;
			B_dev.block(i, j, k)[2] = chi_B_field; 
			*/
			// // end hardcoded B field ////////////////////

		}
		else
		{
			B_dev.block(i, j, k)[0] = 0.0;
			B_dev.block(i, j, k)[1] = 0.0;
			B_dev.block(i, j, k)[2] = 0.0;
		}
		
		
		if (zero_velocities)
		{
			v_b_dev.block(i, j, k)[0] = 0.0;					
			v_b_dev.block(i, j, k)[1] = 0.0;					
			v_b_dev.block(i, j, k)[2] = 0.0;	
		}
		else
		{
			// convert to spherical coordinates
			auto v_spherical = convert_cartesian_to_spherical(tmp_vector[6], 
														      tmp_vector[7], 
														      tmp_vector[8]);		
			v_b_dev.block(i, j, k)[0] = v_spherical[0];					
			v_b_dev.block(i, j, k)[1] = v_spherical[1];					
			v_b_dev.block(i, j, k)[2] = v_spherical[2];	
		}
				
		{
		// continuum 
		Real kappa, sigma, epsilon;
		read_single_node_triple_field(f_cul, i_global, j_global, k_reverse, kappa, sigma, epsilon);
		for (int n = 0; n < N_nu_; ++n)
		{			
			// hardcoded to 0.0 as in PORTA
			// sigma_dev.block(   i, j, k)[n] = 0.0;		
			// k_c_dev.block(     i, j, k)[n] = tmp_vector[12];		
			// eps_c_th_dev.block(i, j, k)[n] = tmp_vector[13];

			sigma_dev.block(   i, j, k)[n] = sigma;
			k_c_dev.block(     i, j, k)[n] = kappa;
			eps_c_th_dev.block(i, j, k)[n] = epsilon;	

			// TODO: read sigma from file .back, da controllare nell'input (solo un valore)
			// read_single_node_single_field(filename_qel,i_global,j_global,k_reverse);				
		}		
		}	
	});
	
	// close files
	MPI_CHECK(MPI_File_close(&fh));		
	MPI_CHECK(MPI_File_close(&f_cul));		
	MPI_CHECK(MPI_File_close(&f_qel));		
}


Real RT_problem::read_single_node_single_field(MPI_File input_file, const int i, const int j, const int k){

	// Output
    Real output;

    // Compute jump
    const int jump = N_x_ * (N_y_ * k + j) + i;

    // Jump to data of interest
    MPI_CHECK(MPI_File_seek(input_file, jump * sizeof(Real), MPI_SEEK_SET));

    // 0 // epsilon
    MPI_CHECK(MPI_File_read(input_file, &output, 1, MPI_DOUBLE, MPI_STATUS_IGNORE)); 
    
    return output;
}

void RT_problem::read_single_node_triple_field(MPI_File input_file, const int i, const int j, const int k, Real &kappa, Real &sigma, Real &epsilon){

	// Output
    Real output;

    // Compute jump
    const int jump = N_x_ * (N_y_ * k + j) + i;

    // Jump to data of interest
    MPI_CHECK(MPI_File_seek(input_file, 3 * jump * sizeof(Real), MPI_SEEK_SET));
    MPI_CHECK(MPI_File_read(input_file, &output, 1, MPI_DOUBLE, MPI_STATUS_IGNORE)); 
    
    kappa = output;

	// Jump to data of interest
	MPI_CHECK(MPI_File_seek(input_file,  (3 * jump + 1) * sizeof(Real) , MPI_SEEK_SET));
	MPI_CHECK(MPI_File_read(input_file, &output, 1, MPI_DOUBLE, MPI_STATUS_IGNORE)); 
	
	sigma = output;

	// Jump to data of interest
	MPI_CHECK(MPI_File_seek(input_file,  (3 * jump + 2) * sizeof(Real), MPI_SEEK_SET));
	MPI_CHECK(MPI_File_read(input_file, &output, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));

	epsilon = output;

}


// read atom and grid quantities 
void RT_problem::read_3D(const char* filename){

	if (mpi_rank_ == 0) std::cout << "Reading PORTA input from " << filename << std::endl;

	const bool zero_velocities = true;
	if (mpi_rank_ == 0 and zero_velocities) std::cout << "ZERO velocities HARDCODED!" << std::endl;

	// reading atom and grids from file
	MPI_File fh;
	MPI_CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh));

	// buffers
	Real entry;
	char c;

    // some irrelevant data
	for (int i = 0; i < PMD_MAIN_HEADER1 - 48; i++) { MPI_CHECK(MPI_File_read_all(fh, &c, 1, MPI_CHAR, MPI_STATUS_IGNORE)); }

	// TODO there are not used 
	double Lx, Ly, Lz, x_origin, y_origin, z_origin; // 48 in previuos line for these
	
	MPI_CHECK(MPI_File_read_all(fh, &Lx, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &Ly, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &Lz, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));

	MPI_CHECK(MPI_File_read_all(fh, &x_origin, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &y_origin, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &z_origin, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));

	MPI_CHECK(MPI_File_read_all(fh, &N_x_, 1, MPI_INT, MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &N_y_, 1, MPI_INT, MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &N_z_, 1, MPI_INT, MPI_STATUS_IGNORE));
	
	// WARNING: valid for Lx = Ly = constant 
	L_ = 1e-5 * Lx/N_x_; // conversion from cm to km

	// some irrelevant data (x,y coordinates)	
	int skip_size = 16384 * sizeof(double);
	MPI_CHECK(MPI_File_seek(fh, skip_size, MPI_SEEK_CUR));

	// get z grid (depth)
	for (int i = 0; i < N_z_; ++i)
	{		
		MPI_CHECK(MPI_File_read_all(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
		 
		entry *= 1e-5; // conversion from cm to km			
		depth_grid_.push_back(entry);      		
	}
	
	// reverse to get correct ordering
	reverse(depth_grid_.begin(), depth_grid_.end());

	// some irrelevant data (extra z coordinates)
	skip_size = (8192 - N_z_) * sizeof(double);	
	MPI_CHECK(MPI_File_seek(fh, skip_size, MPI_SEEK_CUR));

	// inclinations and azimuths per octant 	
	MPI_CHECK(MPI_File_read_all(fh, &N_theta_, 1, MPI_INT , MPI_STATUS_IGNORE));
	MPI_CHECK(MPI_File_read_all(fh, &N_chi_  , 1, MPI_INT , MPI_STATUS_IGNORE));

	// convert from octant to total 
	N_theta_ *= 2;
	N_chi_   *= 4;

	if (mpi_rank_ == 0) std::cout << "N_theta = " << N_theta_ << ", N_chi = " << N_chi_ << ", from PORTA input." << std::endl;
		
	// some irrelevant data	
	skip_size = PMD_MAIN_HEADER2 - 196608 - 16;
	MPI_CHECK(MPI_File_seek(fh, skip_size, MPI_SEEK_CUR));

    int module_head_size;

	// Size of each grid node	
	MPI_CHECK(MPI_File_read_all(fh, &module_head_size, 1, MPI_INT, MPI_STATUS_IGNORE));

	// Size of each grid node	
	MPI_CHECK(MPI_File_read_all(fh, &node_size_, 1, MPI_INT, MPI_STATUS_IGNORE));

    // Jump to data
    header_size_ = PMD_MAIN_HEADER1 + PMD_MAIN_HEADER2 + 12 + module_head_size;

	// some irrelevant data    
    MPI_CHECK(MPI_File_seek(fh, TWOLEVEL_HEADER1, MPI_SEEK_CUR));

	// reading atomic data    
    MPI_CHECK(MPI_File_read_all(fh, &mass_,  1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &Aul_,   1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &Eu_,    1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &Jl2_,   1, MPI_INT,    MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &Ju2_,   1, MPI_INT,    MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &gl_,    1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &gu_,    1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_File_read_all(fh, &T_ref_, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));

    // // change Eu units to [cm-1]     
    // Eu_ /= c_ * h_; 
    
    // hardcoded from FAL-C 
    Eu_ = 23652.304;    
    
  	// some irrelevant data double temp[ny][nx];    matrix of ground (iz=0) for Planckian boundary
	skip_size = (N_x_ * N_y_) * sizeof(double);
	MPI_CHECK(MPI_File_seek(fh, skip_size, MPI_SEEK_CUR));		

	// set angualr grids and sizes and print
	set_theta_chi_grids(N_theta_, N_chi_);
	set_sizes();

	/////////////////// set sizes
	N_nu_ = nu_grid_.size();
	block_size_ = (PetscInt) 4 * N_nu_ * N_theta_ * N_chi_;
	tot_size_   = (PetscInt) N_s_ * block_size_;	
	
	// create space grid
	space_grid_ = std::make_shared<Grid_t>();
	
	// menage grid distribution and init
	set_grid_partition();	
	space_grid_->init(MPI_COMM_WORLD, {N_x_, N_y_, N_z_}, {1, 1, 0},
									 {mpi_size_x_, mpi_size_y_, mpi_size_z_}, use_ghost_layers_); 

	// // TEST
	// space_grid_->init(MPI_COMM_WORLD, {N_x_, N_y_, N_z_}, {1, 1, 0}, {}, use_ghost_layers_); 
		
	// init fields
	allocate_fields();				

	// init atmospheric quantities 
	allocate_atmosphere();	

	// load devices for parallel data structures
	auto T_dev   = T_   ->view_device();
	auto a_dev   = a_   ->view_device();
	auto Nl_dev  = Nl_  ->view_device();
	auto Cul_dev = Cul_ ->view_device();
	auto B_dev   = B_   ->view_device();
	auto v_b_dev = v_b_ ->view_device();
	auto xi_dev  = xi_  ->view_device();
	auto Qel_dev = Qel_ ->view_device();
	auto D2_dev  = D2_  ->view_device();

	auto sigma_dev    = sigma_    ->view_device();
	auto k_c_dev      = k_c_      ->view_device();
	auto eps_c_th_dev = eps_c_th_ ->view_device();
	auto epsilon_dev  = epsilon_  ->view_device();

	auto Doppler_width_dev = Doppler_width_->view_device();

	auto g_dev = space_grid_->view_device();

	// some constants hardcoded	
	if (mpi_rank_ == 0) std::cout << "WARNING: hardcoding quantities for PORTA input read!" << std::endl;
	const double Aul_RH = 2.17674e8;    
	const double mass_real_RH = 6.65511e-23; // hardcoded		

	// hardcoded xi
	std::array<double, 133> xi_vec = { 2.0000000000000000e+00, 2.0000000000000000e+00, 2.0000000000000000e+00, 2.0000000000000000e+00, 2.0000000000000000e+00, 2.0000000000000000e+00, 2.0000000000000000e+00, 1.9886373912499999e+00, 1.9773626719999999e+00, 1.9419696239999999e+00, 1.8600738539999999e+00, 1.7778424079999999e+00, 1.6824751599999999e+00, 1.5795074739999999e+00, 1.4744281220000000e+00, 1.3537171600000000e+00, 1.2327272599999999e+00, 1.1129887319999998e+00, 9.9467995799999975e-01, 8.7645056999999993e-01, 7.8089990399999998e-01, 6.8735903999999992e-01, 6.0959461199999998e-01, 5.6051231999999995e-01, 5.1142826399999997e-01, 4.9363330559999996e-01, 4.8427622399999998e-01, 4.7908703999999996e-01, 4.8846167039999999e-01, 4.9785317760000003e-01, 5.1861966000000015e-01, 5.4516093000000010e-01, 5.8094278171428582e-01, 6.2224760000000012e-01, 6.7522791428571438e-01, 7.3142138000000012e-01, 8.0070458000000022e-01, 8.7089629200000007e-01, 9.5220700000000047e-01, 1.0332738400000001e+00, 1.1213362909090909e+00, 1.2152569890909093e+00, 1.3094004363636367e+00, 1.4048512000000004e+00, 1.5007467520000004e+00, 1.5978336000000002e+00, 1.6978324000000005e+00, 1.7978304000000003e+00, 1.8947160320000005e+00, 1.9907144960000005e+00, 2.0867125760000005e+00, 2.1827118080000001e+00, 2.2787102720000001e+00, 2.3705575253333340e+00, 2.4612234666666670e+00, 2.5518886826666671e+00, 2.6425844611764711e+00, 2.7343480658823536e+00, 2.8261127717647065e+00, 2.9178771105882353e+00, 3.0096418164705878e+00, 3.0936051408695651e+00, 3.1753403478260873e+00, 3.2570791513043482e+00, 3.3388189356521742e+00, 3.4205534886956523e+00, 3.5020644897959188e+00, 3.5755338775510208e+00, 3.6489997387755104e+00, 3.7224691265306125e+00, 3.7959382204081633e+00, 3.8687040000000006e+00, 3.9367072639999998e+00, 4.0047040000000003e+00, 4.0727037280000005e+00, 4.1407040000000004e+00, 4.2084453608247427e+00, 4.2744210474226803e+00, 4.3404041237113411e+00, 4.4063803381443298e+00, 4.4723597195876286e+00, 4.5390595657142860e+00, 4.6076277028571422e+00, 4.6761985828571433e+00, 4.7447738514285716e+00, 4.8133414400000003e+00, 4.8820457066666672e+00, 4.9553796266666668e+00, 5.0287088533333337e+00, 5.1020427733333342e+00, 5.1753755200000002e+00, 5.2502169904761908e+00, 5.3264074666666668e+00, 5.4025979428571436e+00, 5.4787847619047625e+00, 5.5549752380952384e+00, 5.6335312941176472e+00, 5.7182371764705886e+00, 5.8029430588235300e+00, 5.8876489411764705e+00, 5.9725229090909098e+00, 6.0634269090909090e+00, 6.1543359999999998e+00, 6.2481976369230772e+00, 6.3497304123076930e+00, 6.4512745600000008e+00, 6.5531384216216226e+00, 6.6666519351351363e+00, 6.7818653538461549e+00, 6.9049422769230784e+00, 7.0308618105263161e+00, 7.1716484000000005e+00, 7.3142758956521750e+00, 7.4573813333333341e+00, 7.6035485714285720e+00, 7.7321138285714301e+00, 7.8826903272727300e+00, 8.3409513513513822e+00, 1.1900000000000000e+01, 1.1900000000000000e+01, 1.1900000000000000e+01, 1.1900000000000000e+01, 1.1900000000000000e+01, 1.1900000000000000e+01, 1.1900000000000000e+01, 1.1900000000000000e+01, 1.1900000000000000e+01, 1.1900000000000000e+01, 1.1900000000000000e+01, 1.1900000000000000e+01, 1.1900000000000000e+01, 1.1900000000000000e+01, 1.1900000000000000e+01};

	// fill field 
	sgrid::parallel_for("READ-ATM1D", space_grid_->md_range(), SGRID_LAMBDA(int i, int j, int k) {

		// get global indeces form local ones
		const int i_global  = g_dev.start[0] + i - g_dev.margin[0];
		const int j_global  = g_dev.start[1] + j - g_dev.margin[1];
		const int k_global  = g_dev.start[2] + k - g_dev.margin[2];
		
		// reversing z index because of input ordering 
		const int k_reverse = (N_z_ - k_global - 1);  

		auto tmp_vector = read_single_node(fh,i_global,j_global,k_reverse);			

		epsilon_dev.ref(i,j,k) = tmp_vector[0];		
		Cul_dev.ref(i,j,k)     = tmp_vector[1];		
		T_dev.ref(i,j,k)       = tmp_vector[2];			
		Nl_dev.ref(i,j,k)      = tmp_vector[9];		
		a_dev.ref(i,j,k)       = tmp_vector[10];		
		D2_dev.ref(i,j,k)      = tmp_vector[11];
		
		// hardcoding xi
		xi_dev.ref(i,j,k) = 0; // with conversion to cm/s

		// RH Doppler_width to compute Qel
		const double xi = 1e5 * xi_vec[k_reverse];; // with conversion to cm/s
		const double Dw_RH = Eu_ * std::sqrt(xi * xi + 2 * k_B_ * T_dev.ref(i,j,k) / mass_real_RH);		

		// compute Qel
		Qel_dev.ref(i,j,k) = a_dev.ref(i,j,k) * (4 * PI * Dw_RH) - Aul_RH;
		
		// convert to spherical coordinates
		auto B_spherical = convert_cartesian_to_spherical(tmp_vector[3], 
														  tmp_vector[4],
														  tmp_vector[5]);	
		if (use_magnetic_field_)
		{								
			B_dev.block(i, j, k)[0] = B_spherical[0] * 1399600.0; // converting to Larmor frequency					
			B_dev.block(i, j, k)[1] = B_spherical[1]; 					
			B_dev.block(i, j, k)[2] = B_spherical[2]; 

			// // /*  hardcoded B field */ ////////////////////
			/*
			const double B_field_hardcoded = 1000.0; // [Gauss]
			const double theta_B_field = 1.5707963268; // [rad]
			const double chi_B_field = 0.0; // [rad]
			if ( mpi_rank_ == 0 and i == 0 and j == 0 and k == 0) std::cout << "WARNING: HARDCODED B FIELD: of:    " <<  B_field_hardcoded << " Gauss" << std::endl;
			if ( mpi_rank_ == 0 and i == 0 and j == 0 and k == 0) std::cout << "WARNING: HARDCODED B FIELD, theta: " <<  theta_B_field << " rad" << std::endl;
			if ( mpi_rank_ == 0 and i == 0 and j == 0 and k == 0) std::cout << "WARNING: HARDCODED B FIELD, chi:   " <<  chi_B_field << " rad" << std::endl;

			B_dev.block(i, j, k)[0] = GAUSS_TO_LARMOR_FREQUENCY(B_field_hardcoded) ; // converting to Larmor frequency					
			B_dev.block(i, j, k)[1] = theta_B_field;
			B_dev.block(i, j, k)[2] = chi_B_field; 
			*/
			// // end hardcoded B field ////////////////////

		}
		else
		{
			B_dev.block(i, j, k)[0] = 0.0;
			B_dev.block(i, j, k)[1] = 0.0;
			B_dev.block(i, j, k)[2] = 0.0;
		}
		
		
		if (zero_velocities)
		{
			v_b_dev.block(i, j, k)[0] = 0.0;					
			v_b_dev.block(i, j, k)[1] = 0.0;					
			v_b_dev.block(i, j, k)[2] = 0.0;	
		}
		else
		{
			// convert to spherical coordinates
			auto v_spherical = convert_cartesian_to_spherical(tmp_vector[6], 
														      tmp_vector[7], 
														      tmp_vector[8]);		
			v_b_dev.block(i, j, k)[0] = v_spherical[0];					
			v_b_dev.block(i, j, k)[1] = v_spherical[1];					
			v_b_dev.block(i, j, k)[2] = v_spherical[2];	
		}
		
		
		// continuum 
		for (int n = 0; n < N_nu_; ++n)
		{			
			// hardcoded to 0.0 as in PORTA
			sigma_dev.block(   i, j, k)[n] = 0.0;		
			k_c_dev.block(     i, j, k)[n] = tmp_vector[12];		
			eps_c_th_dev.block(i, j, k)[n] = tmp_vector[13];					
		}			
	});

	// close file
	MPI_CHECK(MPI_File_close(&fh));		
}


// read single node from Tanausu
std::vector<Real> RT_problem::read_single_node(MPI_File fh, const int i, const int j, const int k){

    // Output
    std::vector<Real> output;

    // Could be global
	const int L_LIMIT = (Jl2_+1)*(Jl2_+1);
	const int U_LIMIT = (Ju2_+1)*(Ju2_+1);
	const int NJKQ = 9;

	// some temporary variables
	Real entry, atomic_density, rho00l, Nl, Cul;

    // Compute jump
    const int jump = header_size_ + node_size_ * (N_x_ * (N_y_ * k + j) + i);

    // Jump to data of interest
    MPI_CHECK(MPI_File_seek(fh, jump, MPI_SEEK_SET));	

    // 0 // epsilon
    MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE)); 
    output.push_back(entry);		
    
    // 1 // get Cul from epsilon
    Cul = Aul_ * entry / (1.0 - entry); 
    output.push_back(Cul);
    
    // 2 // temperature    
    MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    output.push_back(entry);

   	// atomic density	
    MPI_CHECK(MPI_File_read(fh, &atomic_density, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));	

    // 3 // Bx	    
    MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));	
    output.push_back(entry);  
    // 4 // By    
    MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    output.push_back(entry);		
    // 5 // Bz    
    MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    output.push_back(entry);		
    // 6 // vx	    
    MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));	
    output.push_back(entry);		
    // 7 // vx	    
    MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    output.push_back(entry);		
    // 8 // vz
    MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));
    output.push_back(entry);		

    // density matrix components of the lower level:
    if (L_LIMIT != 1) std::cout << "WARNING: reading 3D input L_LIMIT is not unity!" << std::endl;

    for (int jj = 0; jj < L_LIMIT; jj++) 
    {			        
        MPI_CHECK(MPI_File_read(fh, &rho00l, 1, MPI_DOUBLE, MPI_STATUS_IGNORE)); // real component	
        MPI_CHECK(MPI_File_read(fh, &entry , 1, MPI_DOUBLE, MPI_STATUS_IGNORE)); // im components				
    }    
    
    // recover populations 
    // 9
    Nl = atomic_density * sqrt(Jl2_ + 1) * rho00l;
    output.push_back(Nl);	
   	
    // density matrix components of the upper level:
    for (int jj = 0; jj < U_LIMIT; jj++) 
    {        
        MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE)); // real component	
        MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));	// im components							
    }
    
    // components of the J^K_Q tensor:    
    for (int jj = 0; jj < NJKQ; jj++) MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));	
    

    // rest of the MHD quantities 
    // 10 // Voight parameter a    
	MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));							
    output.push_back(entry);		
    // 11 // delta2 (depolarizing collisional rate)    
    MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));							
    entry *= Aul_; // delta2 = D2/Aul, PORTA uses delta2 here and conversion is needed
    output.push_back(entry);		
    // 12 // continuum opacity    
    MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));							
    output.push_back(entry);		
    // 13 // continuum emissivity
    MPI_CHECK(MPI_File_read(fh, &entry, 1, MPI_DOUBLE, MPI_STATUS_IGNORE));							
    output.push_back(entry);	

    return output;
}


void RT_problem::read_atom(input_string filename){

	if (mpi_rank_ == 0) std::cout << "Reading atom data from " << filename << std::endl;

	std::ifstream myFile(filename);
	std::string line;	

	if (not myFile.good()) std::cerr << "\nERROR: File " << filename << " does not exist!\n" << std::endl;

	int counter = 0;

	while(getline(myFile, line))
	{
		std::istringstream lineStream(line);
		
		std::string line_entry;

		switch(counter) {
    		case 0 : lineStream >> mass_; 
    			     lineStream >> line_entry;
    			     if (line_entry.compare("!Atomic") != 0) std::cout << "WARNING: " << line_entry << " is not !Atomic" << '\n';
            	break;      
    		
    		case 1 : lineStream >> El_;
    				 lineStream >> Jl2_;
    				 lineStream >> gl_;
    				 lineStream >> line_entry;
    				 // if (line_entry.compare("!El[cm-1]") != 0) std::cout << "WARNING: " << line_entry << " is not !El[cm-1]" << '\n';
             	break;
            
            case 2 : lineStream >> Eu_;
    				 lineStream >> Ju2_;
    				 lineStream >> gu_;
    				 lineStream >> line_entry;
    				 // if (line_entry.compare("!Eu[cm-1]") != 0) std::cout << "WARNING: " << line_entry << " is not !Eu[cm-1]" << '\n';
             	break;
            
            case 3 : lineStream >> Aul_;
            		 lineStream >> line_entry;
            		 if (line_entry.compare("!Aul[s-1]") != 0) std::cout << "WARNING: " << line_entry << " is not !Aul[s-1]" << '\n';
             	break;
		}

		counter++;
	}	
}


void RT_problem::read_continumm_1D(input_string filename_sigma, input_string filename_k_c, input_string filename_eps_c_th){

	if (mpi_rank_ == 0) std::cout << "Reading sigma from "    << filename_sigma    << std::endl;
	if (mpi_rank_ == 0) std::cout << "Reading k_c from "      << filename_k_c      << std::endl;
	if (mpi_rank_ == 0) std::cout << "Reading eps_c_th from " << filename_eps_c_th << std::endl;

	std::ifstream myFile_sigma(filename_sigma);
	std::ifstream myFile_k_c(filename_k_c);
	std::ifstream myFile_eps_c_th(filename_eps_c_th);

	std::string line;	

	if (not myFile_sigma.good())    std::cerr << "\nERROR: File " << filename_sigma    << " does not exist!\n" << std::endl;
	if (not myFile_k_c.good())      std::cerr << "\nERROR: File " << filename_k_c      << " does not exist!\n" << std::endl;
	if (not myFile_eps_c_th.good()) std::cerr << "\nERROR: File " << filename_eps_c_th << " does not exist!\n" << std::endl;

	bool first_line = true;

	auto sigma_dev    = sigma_    ->view_device();
	auto k_c_dev      = k_c_      ->view_device();
	auto eps_c_th_dev = eps_c_th_ ->view_device();

	Real entry;	

	std::vector<Real> sigma_vec;
	std::vector<Real> k_c_vec;
	std::vector<Real> eps_c_th_vec;

	while(getline(myFile_sigma, line))
	{	
		std::stringstream iss(line);
	
		if (not first_line) while (iss >> entry) sigma_vec.push_back(entry);					
					
		first_line = false;
	} 

	first_line = true;

	while(getline(myFile_k_c, line))
	{	
		std::stringstream iss(line);
	
		if (not first_line) while (iss >> entry) k_c_vec.push_back(entry);					
					
		first_line = false;
	} 

	first_line = true;
		
	while(getline(myFile_eps_c_th, line))
	{	
		std::stringstream iss(line);
	
		if (not first_line) while (iss >> entry) eps_c_th_vec.push_back(entry);					
					
		first_line = false;
	} 

	const int N_z_N_nu_ = N_z_ * N_nu_;

	// safety check
	if (sigma_vec.size()    != N_z_N_nu_) std::cout << "WARNING: size mismatch in read_continumm_1D()" << std::endl;
	if (k_c_vec.size()      != N_z_N_nu_) std::cout << "WARNING: size mismatch in read_continumm_1D()" << std::endl;
	if (eps_c_th_vec.size() != N_z_N_nu_) std::cout << "WARNING: size mismatch in read_continumm_1D()" << std::endl;
	
	auto g_dev = space_grid_->view_device();

	// fill field
	sgrid::parallel_for("READ SIGMA", space_grid_->md_range(), SGRID_LAMBDA(int i, int j, int k) {

		int k_global;

		for (int n = 0; n < N_nu_; ++n)
		{
			k_global = N_nu_ * (g_dev.start[2] + k - g_dev.margin[2]) + n;

			sigma_dev.block(   i, j, k)[n] = sigma_vec[k_global];		
			k_c_dev.block(     i, j, k)[n] = k_c_vec[k_global];		
			eps_c_th_dev.block(i, j, k)[n] = eps_c_th_vec[k_global];					
		}			
	});	
}


void RT_problem::read_magnetic_field_1D(input_string filename){

	if (mpi_rank_ == 0) std::cout << "Reading magnetic field from " << filename << std::endl;

	std::ifstream myFile(filename);
	std::string line;	

	if (not myFile.good()) std::cerr << "\nERROR: File " << filename << " does not exist!\n" << std::endl;

	bool first_line = true;

	auto B_dev = B_ ->view_device();

	Real entry;	
	std::string entry_label;		

	std::vector<Real> nu_L_vec;
	std::vector<Real> theta_B_vec;
	std::vector<Real> chi_B_vec;

	while(getline(myFile, line))
	{
		std::istringstream lineStream(line);		
		
		if (first_line) // skip first line 
		{
			lineStream >> entry_label;
			if (entry_label.compare("B[G]") != 0)       std::cout << entry_label << " is not B[G]"       << '\n';
			lineStream >> entry_label;
			if (entry_label.compare("theta[rad]") != 0) std::cout << entry_label << " is not theta[rad]" << '\n';
			lineStream >> entry_label;
			if (entry_label.compare("chi[rad]") != 0)   std::cout << entry_label << " is not chi[rad]"   << '\n';			
		}
		else
		{						
			// read data
			lineStream >> entry;				
			entry *= 1399600.0; // convert to Larmor frequency
			nu_L_vec.push_back(entry);
			lineStream >> entry;
			theta_B_vec.push_back(entry);
			lineStream >> entry;
			chi_B_vec.push_back(entry);				
		}		

		first_line = false;
	} 

	// safety check
	if (nu_L_vec.size()    != N_z_) std::cout << "WARNING: size mismatch (" << nu_L_vec.size()    << ") in read_magnetic_field_1D()" << std::endl;
	if (theta_B_vec.size() != N_z_) std::cout << "WARNING: size mismatch (" << theta_B_vec.size() << ") in read_magnetic_field_1D()" << std::endl;
	if (chi_B_vec.size()   != N_z_) std::cout << "WARNING: size mismatch (" << chi_B_vec.size()   << ") in read_magnetic_field_1D()" << std::endl;
	
	auto g_dev = space_grid_->view_device();

	const bool B_etero = false;

	if (mpi_rank_ == 0 and B_etero) std::cout << "\nWARNING: adding perturbation to magnetic field!" << std::endl;

	// fill field
	sgrid::parallel_for("READ B", space_grid_->md_range(), SGRID_LAMBDA(int i, int j, int k) {
				
		const int i_global = g_dev.start[0] + i - g_dev.margin[0];
		const int j_global = g_dev.start[1] + j - g_dev.margin[1];
		const int k_global = g_dev.start[2] + k - g_dev.margin[2];

		if (B_etero)
		{
			const double x = 2.0 * PI * i_global / (N_x_ - 1.0);
			const double y = 2.0 * PI * j_global / (N_y_ - 1.0);		

			const double B_max = 30.0;	

			B_dev.block(i,j,k)[0] = 1399600.0 * B_max * (1 + cos(x) * cos(y)/2.0);	
			B_dev.block(i,j,k)[1] =     PI * abs(cos(x) * cos(y));
			B_dev.block(i,j,k)[2] = 2 * PI * abs(cos(x) * cos(y));
		}
		else
		{
			B_dev.block(i, j, k)[0] =    nu_L_vec[k_global];					
			B_dev.block(i, j, k)[1] = theta_B_vec[k_global];					
			B_dev.block(i, j, k)[2] =   chi_B_vec[k_global];
		}
	});		
}


void RT_problem::read_bulk_velocity_1D(input_string filename){
	
	if (mpi_rank_ == 0) std::cout << "Reading bulk velocities from " << filename << std::endl;

	std::ifstream myFile(filename);
	std::string line;	

	if (not myFile.good()) std::cerr << "\nERROR: File " << filename << " does not exist!\n" << std::endl;

	bool first_line = true;

	auto v_b_dev = v_b_ ->view_device();
	
	Real entry;
	std::string entry_label;		

	std::vector<Real> v_b_vec;
	std::vector<Real> theta_b_vec;
	std::vector<Real> chi_b_vec;
	
	while(getline(myFile, line))
	{
		std::istringstream lineStream(line);
		
		if (first_line) // skip first line 
		{
			lineStream >> entry_label;
			if (entry_label.compare("V[kms-1]") != 0)   std::cout << entry_label << " is not V[kms-1]"   << '\n';
			lineStream >> entry_label;
			if (entry_label.compare("theta[rad]") != 0) std::cout << entry_label << " is not theta[rad]" << '\n';
			lineStream >> entry_label;
			if (entry_label.compare("chi[rad]") != 0)   std::cout << entry_label << " is not chi[rad]"   << '\n';			
		}
		else
		{				
			// read data
			lineStream >> entry;
			v_b_vec.push_back(1e5 * entry); // conversion to cm
			lineStream >> entry;
			theta_b_vec.push_back(entry);
			lineStream >> entry;
			chi_b_vec.push_back(entry);					
		}	

		first_line = false;
	} 

	// safety check
	if (v_b_vec.size()     != N_z_) std::cout << "WARNING: size mismatch in read_bulk_velocity_1D()" << std::endl;
	if (theta_b_vec.size() != N_z_) std::cout << "WARNING: size mismatch in read_bulk_velocity_1D()" << std::endl;
	if (chi_b_vec.size()   != N_z_) std::cout << "WARNING: size mismatch in read_bulk_velocity_1D()" << std::endl;

	auto g_dev = space_grid_->view_device();

	// fill field
	sgrid::parallel_for("READ BULK-VEL", space_grid_->md_range(), SGRID_LAMBDA(int i, int j, int k) {
		
		const int k_global = g_dev.start[2] + k - g_dev.margin[2];

		v_b_dev.block(i, j, k)[0] =     v_b_vec[k_global];					
		v_b_dev.block(i, j, k)[1] = theta_b_vec[k_global];					
		v_b_dev.block(i, j, k)[2] =   chi_b_vec[k_global];															
    });
}


void RT_problem::read_atmosphere_1D(input_string filename){

	if (mpi_rank_ == 0) std::cout << "Reading atmospheric data from " << filename << std::endl;

	std::ifstream myFile(filename);
	std::string line;	

	if (not myFile.good()) std::cerr << "\nERROR: File " << filename << " does not exist!\n" << std::endl;

	bool first_line = true;

	auto T_dev   = T_   ->view_device();
	auto xi_dev  = xi_  ->view_device();
	auto a_dev   = a_   ->view_device();
	auto Nl_dev  = Nl_  ->view_device();
	// auto Nu_dev   = Nu_  ->view_device();
	auto Cul_dev = Cul_ ->view_device();
	auto Qel_dev = Qel_ ->view_device();
	
	Real entry;
	std::string entry_label;
	
	std::vector<Real> T_vec;
	std::vector<Real> xi_vec;
	std::vector<Real> a_vec;
	std::vector<Real> Nl_vec;
	// std::vector<Real> Nu_vec;
	std::vector<Real> Cul_vec;
	std::vector<Real> Qel_vec;
	
	while(getline(myFile, line))
	{	
		std::istringstream lineStream(line);	

		if (first_line) // skip first line 
		{						
			lineStream >> entry_label;
			if (entry_label.compare("Height[km]") != 0) std::cout << entry_label << " is not Height[km]" << '\n';			
			lineStream >> entry_label;
			if (entry_label.compare("Temp[K]") != 0)    std::cout << entry_label << " is not Temp[K]" << '\n';			
			lineStream >> entry_label;
			if (entry_label.compare("Vmic[cm/s]") != 0) std::cout << entry_label << " is not Vmic[cm/s]" << '\n';		
			lineStream >> entry_label;
			if (entry_label.compare("Damp") != 0)       std::cout << entry_label << " Damp" << '\n';			
			lineStream >> entry_label;
			if (entry_label.compare("Nl[cm-3]") != 0)   std::cout << entry_label << " is not Nl[cm-3]" << '\n';			
			lineStream >> entry_label;
			if (entry_label.compare("Nu[cm-3]") != 0)   std::cout << entry_label << " is not Nu[cm-3]" << '\n';			
			lineStream >> entry_label;
			if (entry_label.compare("Cul[s-1]") != 0)   std::cout << entry_label << " is not Cul[s-1]" << '\n';						
			lineStream >> entry_label;
			if (entry_label.compare("Qel[s-1]") != 0)   std::cout << entry_label << " is not Qel[s-1]" << '\n';			
		}
		else
		{					
			// read data
			lineStream >> entry;
			// depth_grid_.push_back(entry);					
			lineStream >> entry;
			T_vec.push_back(entry);
			lineStream >> entry;
			xi_vec.push_back(entry);
			lineStream >> entry;
			a_vec.push_back(entry);
			lineStream >> entry;
			Nl_vec.push_back(entry);
			lineStream >> entry;
			// Nu_vec.push_back(entry);
			lineStream >> entry;
			Cul_vec.push_back(entry);
			lineStream >> entry;
			Qel_vec.push_back(entry);		
		}		

		first_line = false;
	} 	

	// safety check
	if (T_vec.size()   != N_z_) std::cout << "WARNING: size mismatch in T_vec"  << std::endl;
	if (xi_vec.size()  != N_z_) std::cout << "WARNING: size mismatch in xi_vec" << std::endl;
	if (a_vec.size()   != N_z_) std::cout << "WARNING: size mismatch in a"      << std::endl;
	if (Nl_vec.size()  != N_z_) std::cout << "WARNING: size mismatch in Nl"     << std::endl;
	if (Cul_vec.size() != N_z_) std::cout << "WARNING: size mismatch in Cul"    << std::endl;
	if (Qel_vec.size() != N_z_) std::cout << "WARNING: size mismatch in Qel"    << std::endl;

	auto g_dev = space_grid_->view_device();

	const bool test_etero = false;

	if (mpi_rank_ == 0 and test_etero) std::cout << "WARNING: adding temperature perturbation!" << std::endl;

	// fill field 
	sgrid::parallel_for("READ-ATM1D", space_grid_->md_range(), SGRID_LAMBDA(int i, int j, int k) {

		const int k_global = g_dev.start[2] + k - g_dev.margin[2];

		if (N_x_ > 1 and test_etero)
		{
			const int i_global = g_dev.start[0] + i - g_dev.margin[0];
			const int j_global = g_dev.start[1] + j - g_dev.margin[1];

			const double x = i_global / (N_x_ - 1.0);
			const double y = j_global / (N_y_ - 1.0);

			const double T_ijk = T_vec[k_global] * cos(2.0 * PI * x) * cos(2.0 * PI * y) / 2.0;

			T_dev.ref(i,j,k) = T_vec[k_global] + T_ijk;
		}
		else
		{
			T_dev.ref(i,j,k) = T_vec[k_global];
		}		

		xi_dev.ref( i, j, k) =  xi_vec[k_global];		
		a_dev.ref(  i, j, k) =   a_vec[k_global];		
		Nl_dev.ref( i, j, k) =  Nl_vec[k_global];		
		Cul_dev.ref(i, j, k) = Cul_vec[k_global];		
		Qel_dev.ref(i, j, k) = Qel_vec[k_global];						
	});

	// T_->write("T32.raw");          
}


void RT_problem::read_depth(input_string filename){
	
	if (mpi_rank_ == 0) std::cout << "Reading depth data from " << filename << std::endl;

	std::ifstream myFile(filename);
	std::string line;	

	if (not myFile.good()) std::cerr << "\nERROR: File " << filename << " does not exist!\n" << std::endl;

	bool first_line = true;

	while(getline(myFile, line))
	{
		std::istringstream lineStream(line);
		Real entry;
		std::string entry_label;

		if (first_line) // skip first line 
		{
			lineStream >> entry_label;
			if (entry_label.compare("Height[km]") != 0) std::cout << entry_label << " is not Height[km]" << '\n';						
		}
		else
		{			
			lineStream >> entry;					
			depth_grid_.push_back(entry);
		}		

		first_line = false;
	} 		
}


void RT_problem::read_frequency(input_string filename, const bool use_wavelength){

	if (mpi_rank_ == 0)
	{
		if (use_wavelength)
		{
			std::cout << "Reading frequencies in [A] from " << filename << std::endl;
		}
		else
		{
			std::cout << "Reading frequencies in [s-1] from " << filename << std::endl;
		}
	} 

	std::ifstream myFile(filename);
	std::string line;	

	if (not myFile.good()) std::cerr << "\nERROR: File " << filename << " does not exist!\n" << std::endl;

	bool first_line = true;

	while(getline(myFile, line))
	{
		std::istringstream lineStream(line);
		double entry;
		std::string entry_label;

		if (not first_line) // skip first line 
		{		
			if (use_wavelength)
			{
				lineStream >> entry;	

				// TODO conversion from air to vacuum

				// convert to [s-1]
				entry = 1e8 * c_ / entry;

				nu_grid_.push_back(entry);		
			}
			else
			{
				lineStream >> entry;			
				lineStream >> entry;
				// if(this->mpi_rank_ == 0) {
				// 	 std::cout << "Frequency: " << std::scientific << std::setprecision(15) << entry << std::endl;
				// }
				nu_grid_.push_back(entry);					
			}

			if (entry < 1.0e14) std::cerr << "\nWARNING: frequency: " << entry << " probably not in Herz[s-1]!\n" << std::endl;
		}		

		first_line = false;
	} 
}


void RT_problem::set_sizes(){

	// set disc. parameters from input		
	N_z_     = depth_grid_.size();
	N_nu_	 = nu_grid_.size();
	N_theta_ = theta_grid_.size();
	N_chi_	 = chi_grid_.size();
	N_dirs_  = N_theta_ * N_chi_; 

	N_s_ = N_x_ * N_y_ * N_z_;

	block_size_ = (PetscInt) 4 * N_nu_ * N_theta_ * N_chi_;
	tot_size_   = (PetscInt) N_s_ * block_size_;	

	if (mpi_rank_ == 0 and mpi_size_ > N_s_) std::cerr << "\n========= WARNING: mpi_size > N_s! =========\n" << std::endl;

	if (mpi_size_ > block_size_/4)
	{
		if (mpi_rank_ == 0) std::cout << "\nUsing mpi_size > block_size/4" << std::endl;
	}
	else if ((block_size_ / mpi_size_) % 4 != 0)
	{
		// if (mpi_rank_ == 0) 
		if (mpi_rank_ == 0) 
		{
			std::cerr << "\n========= ERROR: block_size_ / mpi_size_ should be a divisible by 4! =========\n" << std::endl;
			std::cout << "block_size = " << block_size_ << std::endl;
			std::cout << "mpi_size_  = " << mpi_size_   << std::endl;
		}

		throw;
	}	
}

void const RT_problem::print_info(){
	
	if (mpi_rank_ == 0) 		
	{		
		std::cout << "\nmpi_size_x = " << mpi_size_x_ << std::endl;
		std::cout <<   "mpi_size_y = " << mpi_size_y_ << std::endl;
		std::cout <<   "mpi_size_z = " << mpi_size_z_ << std::endl << std::endl;

		std::cout << "\n=========== 2-levels atom parameters ===========\n" << std::endl;		
		std::cout << "Mass = " << mass_ << std::endl;
		std::cout << "El = "   << El_   << std::endl;
		std::cout << "Eu = "   << Eu_   << std::endl;		
		std::cout << "2Jl = "  << Jl2_   << std::endl;
		std::cout << "2Ju = "  << Ju2_   << std::endl;
		std::cout << "gl = "   << gl_   << std::endl;		
		std::cout << "gu = "   << gu_   << std::endl;
		std::cout << "Aul = "  << Aul_  << std::endl;
		
		std::cout << "\n=========== Grids parameters ===========\n" << std::endl;	
		std::cout << "N_x = "     << N_x_     << std::endl;	
		std::cout << "N_y = "     << N_y_     << std::endl;	
		std::cout << "N_z = "     << N_z_     << std::endl;	
		std::cout << "N_s = "     << N_s_     << std::endl;	
		std::cout << "N_theta = " << N_theta_ << std::endl;	
		std::cout << "N_chi = "   << N_chi_   << std::endl;	
		std::cout << "N_dirs = "  << N_dirs_  << std::endl;	
		std::cout << "N_nu = "    << N_nu_    << std::endl;	
		std::cout << "L = "       << L_       << std::endl;			
		
		std::cout << "\ntotal size = " << tot_size_   << std::endl;			
		std::cout << "block size = "   << block_size_ << std::endl;		

		std::cout << "\ntheta grid = [ ";

		for (int i = 0; i < N_theta_; ++i) std::cout << theta_grid_[i] << " ";

		std::cout << "]\nmu grid    = [ ";

		for (int i = 0; i < N_theta_; ++i) std::cout   << mu_grid_[i] << " ";

		std::cout << "]\nchi grid   = [ ";

		for (int i = 0; i < N_chi_; ++i) std::cout   << chi_grid_[i] << " ";

		std::cout << "] " << std::endl;			
	}
}


void RT_problem::allocate_fields(){

		// create fields 
		I_field_ = std::make_shared<Field_t>("I", space_grid_, block_size_); 
		S_field_ = std::make_shared<Field_t>("S", space_grid_, block_size_);

		eta_field_ = std::make_shared<Field_t>("eta", space_grid_, block_size_);
		rho_field_ = std::make_shared<Field_t>("rho", space_grid_, block_size_);

		// necessary?
		I_field_->  allocate_on_device(); 
		S_field_->  allocate_on_device(); 
		eta_field_->allocate_on_device(); 
		rho_field_->allocate_on_device(); 		

		///////////////////////

		if (mpi_rank_ == 0) std::cout << "\nCreating PETSc vector..." << std::endl;		

		PetscErrorCode ierr; 
	
		auto g_dev = space_grid_->view_device();

		local_size_ = block_size_ * g_dev.dim[0] * g_dev.dim[1] * g_dev.dim[2];
		
		ierr = VecCreate(PETSC_COMM_WORLD, &I_vec_);CHKERRV(ierr);	
		ierr = VecSetSizes(I_vec_, local_size_, tot_size_);CHKERRV(ierr);			
		ierr = VecSetFromOptions(I_vec_);CHKERRV(ierr);		
}


// allocate fields in a single direction Omega
void RT_problem::allocate_fields_Omega(){		

	const PetscInt block_size_Omega = 4 * N_nu_;

	// create fields 
	I_field_Omega_ = std::make_shared<Field_t>("I_Omega", space_grid_, block_size_Omega); 
	S_field_Omega_ = std::make_shared<Field_t>("S_Omega", space_grid_, block_size_Omega);

	eta_field_Omega_ = std::make_shared<Field_t>("eta_Omega", space_grid_, block_size_Omega);
	rho_field_Omega_ = std::make_shared<Field_t>("rho_Omega", space_grid_, block_size_Omega);

	// necessary?
	I_field_Omega_->  allocate_on_device(); 
	S_field_Omega_->  allocate_on_device(); 
	eta_field_Omega_->allocate_on_device(); 
	rho_field_Omega_->allocate_on_device(); 		
}


void RT_problem::allocate_atmosphere(){
	
	// create atmospheric quantities 
	D1_  = std::make_shared<Field_t>("D1",   space_grid_);
	D2_  = std::make_shared<Field_t>("D2",   space_grid_);
	Nl_  = std::make_shared<Field_t>("Nl",   space_grid_);
	// Nu_   = std::make_shared<Field_t>("Nu",   space_grid);
	T_   = std::make_shared<Field_t>("T",    space_grid_);
	xi_  = std::make_shared<Field_t>("xi",   space_grid_);
	Cul_ = std::make_shared<Field_t>("Cul",  space_grid_);
	Qel_ = std::make_shared<Field_t>("Qel",  space_grid_);
	a_   = std::make_shared<Field_t>("a",    space_grid_);
	W_T_ = std::make_shared<Field_t>("W_T",  space_grid_);

	// magnetic field 
	B_ = std::make_shared<Field_t>("B", space_grid_, 3); 	

	// bulk velocities, in polar coordinates
	v_b_ = std::make_shared<Field_t>("v_b", space_grid_, 3);    
	
	// quantities depending on position that can be precomputed
	Doppler_width_ = std::make_shared<Field_t>("Doppler_width", space_grid_);
	k_L_           = std::make_shared<Field_t>("k_L", 		    space_grid_);
	epsilon_       = std::make_shared<Field_t>("epsilon_", 		space_grid_);

	// input quantities depending on position and frequency 
	u_        = std::make_shared<Field_t>("u",        space_grid_, N_nu_);  
	sigma_    = std::make_shared<Field_t>("sigma",    space_grid_, N_nu_);
	k_c_      = std::make_shared<Field_t>("k_c",      space_grid_, N_nu_);
	eps_c_th_ = std::make_shared<Field_t>("eps_c_th", space_grid_, N_nu_);
	
	// allocate
	D1_  -> allocate_on_device(); 
	D2_  -> allocate_on_device(); 
	Nl_  -> allocate_on_device(); 
	// Nu_  -> allocate_on_device(); 
	T_   -> allocate_on_device(); 
	xi_  -> allocate_on_device(); 	
	Cul_ -> allocate_on_device(); 
	Qel_ -> allocate_on_device(); 
	a_   -> allocate_on_device(); 
	W_T_ -> allocate_on_device(); 	
	B_   -> allocate_on_device(); 
	v_b_ -> allocate_on_device(); 
		
	Doppler_width_ -> allocate_on_device(); 
	k_L_           -> allocate_on_device(); 
	epsilon_       -> allocate_on_device(); 

	u_        -> allocate_on_device(); 
	sigma_ 	  -> allocate_on_device(); 
	k_c_      -> allocate_on_device(); 
	eps_c_th_ -> allocate_on_device(); 
}


void RT_problem::init_field(Field_ptr_t input_field, const Real input_value){

	auto field_dev = input_field->view_device();

    sgrid::parallel_for("INIT I", space_grid_->md_range(), KOKKOS_LAMBDA(int i, int j, int k) 
    {         
        auto *block = field_dev.block(i, j, k);
         
        for (int b = 0; b < block_size_; ++b) 
        {
        	block[b] = input_value;        	
        }
    });
}


void RT_problem::set_theta_chi_grids(const int N_theta, const int N_chi, const bool double_GL){

	if (mpi_rank_ == 0 and N_theta % 2 != 0) std::cerr << "\n========= WARNING: N_theta odd! =========\n" << std::endl;

	// init theta and mu grids and weights 
    // legendre_rule(N_theta,  0.0, PI, theta_grid_, w_theta_);

    if (double_GL)    	
    {
    	std::vector<Real> mu_1;
    	std::vector<Real> mu_2;

    	std::vector<Real> w_1;
    	std::vector<Real> w_2;

    	legendre_rule(N_theta / 2,  -1.0, 0.0, mu_1, w_1);
    	legendre_rule(N_theta / 2,   0.0, 1.0, mu_2, w_2);

    	// first half
    	for (int i = 0; i < N_theta / 2; ++i)
	    {
	    	mu_grid_.push_back(mu_1[i]);
	    	w_theta_.push_back(w_1[i]);

			theta_grid_.push_back(std::acos(mu_1[i]));		    	
	    }

	    // second half
    	for (int i = 0; i < N_theta / 2; ++i)
	    {
	    	mu_grid_.push_back(mu_2[i]);
	    	w_theta_.push_back(w_2[i]);

			theta_grid_.push_back(std::acos(mu_2[i]));		    	
	    }
    }
    else
    {
    	legendre_rule(N_theta,  -1.0, 1.0, mu_grid_, w_theta_);

	    for (int i = 0; i < N_theta; ++i)
	    {
	    	theta_grid_.push_back(std::acos(mu_grid_[i]));
	    }
    }

    // init equidistant chi grid in [0, 2pi] and trap weights

    if (N_chi % 2 != 0) std::cout << "WARNING: chi grid is not even!" << std::endl;
    if (N_chi % 4 != 0) std::cout << "WARNING: chi grid is not the same for each quadrant!" << std::endl;
    
    const Real delta_chi = 2.0 * PI / N_chi;    

    for (int i = 0; i < N_chi; ++i)
    {
    	// chi_grid_.push_back(i * delta_chi + 0.00001);  // grid axes included
    	chi_grid_.push_back((i + 0.5) * delta_chi);	 // avoiding grid axes   	
    	w_chi_.push_back(delta_chi);    	
    }    
}


std::vector<std::complex<Real> > RT_problem::compute_T_KQ(const int stokes_i, const Real theta, const Real chi){

	// using standard ordering from tables: T_KQ = [T_00 T_10 T_11 T_20 T_21 T_22] 
	// for Q negative T_K-Q = (-1)^Q * bar(T_KQ)

	// output
	std::vector<std::complex<Real> > T_KQ(6, 0.0);

	std::complex<Real> eix  = std::polar(1.0, chi);        // exp(i * chi)
	std::complex<Real> e2ix = std::polar(1.0, 2.0 * chi);  // exp(2 * i * chi)

	// various coeffs
 	std::complex<Real> fa, fb;

	Real c2g = std::cos(2.0 * gamma_);
    Real ct  = std::cos(theta);
    Real ct2 = std::cos(theta) * std::cos(theta);
    Real st  = std::sin(theta);    
    Real s2g = std::sin(2.0 * gamma_);
    Real st2 = std::sin(theta) * std::sin(theta);

    Real k  = 1.0 / (2.0 * std::sqrt(2.0));
    Real k2 = 0.5 * std::sqrt(3.0);

	if (stokes_i == 0)
	{
		T_KQ[0] = std::real(1.0); 
		
 		T_KQ[3] = std::real(k * (3.0 * ct2 - 1.0));
		T_KQ[4] = - k2 * st * ct * eix;
		T_KQ[5] =  0.5 * k2 * st2 * e2ix;
	}
	else if (stokes_i == 1)
	{
		T_KQ[3] = std::real(- 3.0 * k * c2g * st2);

	  	fa.real(c2g * ct);
    	fa.imag(s2g);
    	fb = st * eix;	

		T_KQ[4] = - k2 * fa * fb;

		fa.real(c2g * (1.0 + ct2));
    	fa.imag(2.0 * s2g * ct);

		T_KQ[5] = - 0.5 * k2 * fa * e2ix;
	}
	else if (stokes_i == 2)
	{
		T_KQ[3] = std::real(3.0 * k * s2g * st2);

		fa.real(s2g * ct);
    	fa.imag(-c2g);

		T_KQ[4] = k2 * fa * st * eix;

		fa.real(s2g * (1.0 + ct2));
    	fa.imag(-2.0 * c2g * ct);

		T_KQ[5] = 0.5 * k2 * fa * e2ix;
	}
	else if (stokes_i == 3)
	{
		T_KQ[1] = std::real(sqrt(3.0 / 2.0) * ct);
		T_KQ[2] = - k2 * st * eix;
	}
	else
	{
		std::cerr << "\nERROR: stokes_i is too large" << std::endl; 			
	}
		
	return T_KQ;
}


std::complex<Real> RT_problem::get_TKQi(const int i_stokes, const int K, const int Q, const int j, const int k){

	// checks 	
	if ( i_stokes >= 4) std::cerr  << "\nERROR: i_stokes is too large" << std::endl; 					
	if (K > 2) std::cerr           << "\nERROR: K is too large"        << std::endl; 			
	if (std::abs(Q) > K) std::cerr << "\nERROR: Q is too large"        << std::endl; 	
	if (j >= N_theta_) std::cerr   << "\nERROR: N_theta_ is too large" << std::endl; 	
	if (k >= N_chi_) std::cerr     << "\nERROR: N_chi_ is too large"   << std::endl; 	

	int index = i_stokes * N_theta_ * N_chi_ + j * N_chi_ + k;

	std::complex<Real> T_KQ;

	if (K == 0)
	{
		T_KQ = T_KQ_[index][0];
	}
	else if (K == 1)
	{
		if (Q == 0)
		{
			T_KQ = T_KQ_[index][1];
		}
		else if (Q == 1)
		{
			T_KQ = T_KQ_[index][2];
		}
		else if (Q == -1)
		{
			T_KQ = - 1.0 * std::conj(T_KQ_[index][2]);
		}				
		else { std::cerr << "\nERROR: wrong Q input" << std::endl; }		
	}
	else if (K == 2)
	{
		if (Q == 0)
		{
			T_KQ = T_KQ_[index][3];
		}
		else if (Q == 1)
		{
			T_KQ = T_KQ_[index][4];
		}
		else if (Q == -1)
		{
			T_KQ = - 1.0 * std::conj(T_KQ_[index][4]);
		}
		else if (Q == 2)
		{
			T_KQ = T_KQ_[index][5];
		}
		else if (Q == -2)
		{
			T_KQ = std::conj(T_KQ_[index][5]);
		}
		else { std::cerr << "\nERROR: wrong Q input" << std::endl; }		
	}
	else
	{
		std::cerr << "\nERROR: wrong K input" << std::endl; 
	}

	return T_KQ;
}


// for data structure with a single direction
std::complex<Real> RT_problem::get_TKQi(const std::vector<std::complex<Real>> T_KQ_i, const int K, const int Q)
{
	// checks 		
	if (K > 2) std::cerr           << "\nERROR: K is too large"        << std::endl; 			
	if (std::abs(Q) > K) std::cerr << "\nERROR: Q is too large"        << std::endl; 		

	std::complex<Real> T_KQ;

	if (K == 0)
	{
		T_KQ = T_KQ_i[0];
	}
	else if (K == 1)
	{
		if (Q == 0)
		{
			T_KQ = T_KQ_i[1];
		}
		else if (Q == 1)
		{
			T_KQ = T_KQ_i[2];
		}
		else if (Q == -1)
		{
			T_KQ = - 1.0 * std::conj(T_KQ_i[2]);
		}				
		else { std::cerr << "\nERROR: wrong Q input" << std::endl; }		
	}
	else if (K == 2)
	{
		if (Q == 0)
		{
			T_KQ = T_KQ_i[3];
		}
		else if (Q == 1)
		{
			T_KQ = T_KQ_i[4];
		}
		else if (Q == -1)
		{
			T_KQ = - 1.0 * std::conj(T_KQ_i[4]);
		}
		else if (Q == 2)
		{
			T_KQ = T_KQ_i[5];
		}
		else if (Q == -2)
		{
			T_KQ = std::conj(T_KQ_i[5]);
		}
		else { std::cerr << "\nERROR: wrong Q input" << std::endl; }		
	}
	else
	{
		std::cerr << "\nERROR: wrong K input" << std::endl; 
	}

	return T_KQ;
}


void RT_problem::set_eta_and_rhos(){

	auto eta_dev = eta_field_->view_device();
	auto rho_dev = rho_field_->view_device();

	auto a_dev   = a_    ->view_device();
	auto u_dev   = u_    ->view_device();
	auto k_L_dev = k_L_  ->view_device();
	auto k_c_dev = k_c_  ->view_device();	
	auto B_dev   = B_    ->view_device();	
	auto v_b_dev = v_b_  ->view_device();
	
	auto Doppler_width_dev = Doppler_width_->view_device();

    sgrid::parallel_for("INIT ETA-RHO", space_grid_->md_range(), SGRID_LAMBDA(int i, int j, int k) 
    {         
        auto *block_eta = eta_dev.block(i, j, k);
        auto *block_rho = rho_dev.block(i, j, k);
        
        auto *u   =   u_dev.block(i, j, k);		
		auto *k_c = k_c_dev.block(i, j, k);		
		auto *B   =   B_dev.block(i, j, k);       
        auto *v_b = v_b_dev.block(i, j, k);                        
                        
        // assign some variables for readability
        Real theta_v_b = v_b[1];
        Real chi_v_b   = v_b[2];

        Real nu_L    = B[0];
        Real theta_B = B[1];
        Real chi_B   = B[2];

        Real Doppler_width = Doppler_width_dev.ref(i,j,k);
    	Real k_L           = k_L_dev.ref(i,j,k);
		Real a             = a_dev.ref(i,j,k);

		// init rotation matrix
        Rotation_matrix R(0.0, -theta_B, -chi_B);
        
        // indeces
        std::vector<int> local_idx;
        int j_theta, k_chi, n_nu;
         
        for (int b = 0; b < block_size_; b = b + 4) 
        {        	        	
        	block_rho[b + 1] = 0;

			local_idx = block_to_local(b);

        	j_theta = local_idx[0];
        	k_chi   = local_idx[1];
        	n_nu    = local_idx[2];

        	const Real theta = theta_grid_[j_theta];
			const Real chi   = chi_grid_[k_chi];	

			const Real coeff  =  k_L / (std::sqrt(PI) * Doppler_width);
			const Real coeff2 = nu_L / Doppler_width; 

			const std::complex<Real> a_damp(0.0, a);

			// for reduced frequency
			const Real v_dot_Omega = v_b[0] * ( cos(theta_v_b) * cos(theta) + sin(theta_v_b) * sin(theta) * cos(chi - chi_v_b));
			const Real u_b = nu_0_ * v_dot_Omega / (c_ * Doppler_width);

			const Real u_red = u[n_nu] + u_b;			
			
			for (int K = 0; K < 3; ++K)
			{
				const double coeff_K = coeff * std::sqrt(3.0 * (2.0 * double(K) + 1.0));
	
				for (int Mu2 = - Ju2_; Mu2 < Ju2_ + 1; Mu2 += 2)
				{
					for (int Ml2 = - Jl2_; Ml2 < Jl2_ + 1; Ml2 += 2)
					{										
						if (std::abs(Mu2 - Ml2) <= 2) 
						{
							const int q2 = Ml2 - Mu2;

				      		const double W3J1 = W3JS(Ju2_, Jl2_, 2,-Mu2, Ml2, -q2);  
				      		const double W3J2 = W3JS(2, 2, 2 * K, q2, -q2, 0); 

							const double fact = coeff_K * std::pow(-1.0, double(q2) / 2.0 + 1.0) * std::pow(W3J1, 2) * W3J2;								   
							
		      				const double um = coeff2 * (gu_ * (double(Mu2) / 2.0) - gl_ * (double(Ml2) / 2.0)) + u_red;
					
							for (int Q = -K; Q < K + 1; ++Q)
							{			
								const std::complex<double> faddeva = Faddeeva::w(um + a_damp);
				        		const auto D_KQQ                   = std::conj(R(K, 0, Q));

				        		const auto fact_re = fact * std::real(faddeva) * D_KQQ;
				        		const auto fact_im = fact * std::imag(faddeva) * D_KQQ;		
				        		
								// etas							
								block_eta[b    ] += std::real(fact_re * get_TKQi(0, K, Q, j_theta, k_chi)); 
								block_eta[b + 1] += std::real(fact_re * get_TKQi(1, K, Q, j_theta, k_chi));
								block_eta[b + 2] += std::real(fact_re * get_TKQi(2, K, Q, j_theta, k_chi));
								block_eta[b + 3] += std::real(fact_re * get_TKQi(3, K, Q, j_theta, k_chi));														
															
								// rhos
								block_rho[b + 1] += std::real(fact_im * get_TKQi(1, K, Q, j_theta, k_chi));
								block_rho[b + 2] += std::real(fact_im * get_TKQi(2, K, Q, j_theta, k_chi));
								block_rho[b + 3] += std::real(fact_im * get_TKQi(3, K, Q, j_theta, k_chi));									

							}
						}
					}		
				}
        	}

        	if (enable_continuum_) block_eta[b] += k_c[n_nu];       

        	// if (i == 0 and j == 0 and g_dev.global_coord(2, k) == 0 and b >= block_size_ - 4 * N_nu_) 
        	// {
        	// 	// std::cout << "k = "     <<   k      << std::endl; 
        	// 	std::cout <<   block_eta[b]       << std::endl; 
        	// }	
			// std::cout << "block_eta[b + 1] = " <<   block_eta[b + 1]   << std::endl; 
			// std::cout << "block_eta[b + 2] = " <<   block_eta[b + 2]   << std::endl; 
			// std::cout << "block_eta[b + 3] = " <<   block_eta[b + 3]   << std::endl; 
        	// std::cout << "block_rho[b + 1] = " <<   block_rho[b + 1]   << std::endl; 
        	// std::cout << "block_rho[b + 2] = " <<   block_rho[b + 2]   << std::endl; 
        	// std::cout << "block_rho[b + 3] = " <<   block_rho[b + 3]   << std::endl; 

        	// sanity checks
        	if (block_eta[b] == 0) std::cerr << "\nWARNING: zero eta_I!"     << std::endl; 
        	if (block_eta[b] < 0)  std::cerr << "\nWARNING: negative eta_I!" << std::endl; 		
        	if (block_rho[b] < 0)  std::cerr << "\nWARNING: negative rho_I!" << std::endl; 	      

        	if (isnan(block_eta[b    ])) std::cerr << "\nWARNING: eta_I = NaN!" << std::endl; 
        	if (isnan(block_eta[b + 1])) std::cerr << "\nWARNING: eta_Q = NaN!" << std::endl; 
        	if (isnan(block_eta[b + 2])) std::cerr << "\nWARNING: eta_U = NaN!" << std::endl; 
        	if (isnan(block_eta[b + 3])) std::cerr << "\nWARNING: eta_V = NaN!" << std::endl;         	
        	if (isnan(block_rho[b + 1])) std::cerr << "\nWARNING: rho_Q = NaN!" << std::endl; 
        	if (isnan(block_rho[b + 2])) std::cerr << "\nWARNING: rho_U = NaN!" << std::endl; 
        	if (isnan(block_rho[b + 3])) std::cerr << "\nWARNING: rho_V = NaN!" << std::endl;         	      	        	
        } 	
    });	

	// debug			
	// const Real dichroism_module = std::sqrt(etas_and_rhos[1] * etas_and_rhos[1] + etas_and_rhos[2] * etas_and_rhos[2] + etas_and_rhos[3] * etas_and_rhos[3]);	
	// if (etas_and_rhos[0] < dichroism_module) dichroism_warning = true;				
	// if (dichroism_warning) std::cout << "\nWARNING: eta_I < eta! (Eq. (7) Gioele Paganini 2018, Part III)" << std::endl; 		
}


void RT_problem::set_eta_and_rhos_Omega(const Real theta, const Real chi){	

	// vector with KQ components for each stokes profile
	std::vector< std::vector<std::complex<Real> > > T_KQ(4);

	for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
	{
		T_KQ[i_stokes] = compute_T_KQ(i_stokes, theta, chi);		
	}	

	auto eta_dev = eta_field_Omega_->view_device();
	auto rho_dev = rho_field_Omega_->view_device();

	auto a_dev   = a_    ->view_device();
	auto u_dev   = u_    ->view_device();
	auto k_L_dev = k_L_  ->view_device();
	auto k_c_dev = k_c_  ->view_device();	
	auto B_dev   = B_    ->view_device();	
	auto v_b_dev = v_b_  ->view_device();
	
	auto Doppler_width_dev = Doppler_width_->view_device();	
	
    sgrid::parallel_for("INIT ETA-RHO", space_grid_->md_range(), SGRID_LAMBDA(int i, int j, int k) 
    {             	
        auto *block_eta = eta_dev.block(i, j, k);
        auto *block_rho = rho_dev.block(i, j, k);
        
        auto *u   =   u_dev.block(i, j, k);		
		auto *k_c = k_c_dev.block(i, j, k);		
		auto *B   =   B_dev.block(i, j, k);       
        auto *v_b = v_b_dev.block(i, j, k);           
                        
        // assign some variables for readability
        Real theta_v_b = v_b[1];
        Real chi_v_b   = v_b[2];

        Real nu_L    = B[0];
        Real theta_B = B[1];
        Real chi_B   = B[2];

        Real Doppler_width = Doppler_width_dev.ref(i,j,k);
    	Real k_L           = k_L_dev.ref(i,j,k);
		Real a             = a_dev.ref(i,j,k);

		// init rotation matrix
        Rotation_matrix R(0.0, -theta_B, -chi_B);
        
        // indeces
        int b;
        
        for (int n_nu = 0; n_nu < N_nu_; n_nu++) 
        {        	        	        				
			// index
			b = 4 * n_nu;		

			// init to zero
			for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
			{
				block_eta[b + i_stokes] = 0;				
				block_rho[b + i_stokes] = 0;								
			}

			const Real coeff  =  k_L / (std::sqrt(PI) * Doppler_width);
			const Real coeff2 = nu_L / Doppler_width; 

			const std::complex<Real> a_damp(0.0, a);

			// for reduced frequency
			const Real v_dot_Omega = v_b[0] * ( cos(theta_v_b) * cos(theta) + sin(theta_v_b) * sin(theta) * cos(chi - chi_v_b));
			const Real u_b = nu_0_ * v_dot_Omega / (c_ * Doppler_width);

			const Real u_red = u[n_nu] + u_b;

			for (int K = 0; K < 3; ++K)
			{				
				const double coeff_K = coeff * std::sqrt(3.0 * (2.0 * double(K) + 1.0));
					
				for (int Mu2 = - Ju2_; Mu2 < Ju2_ + 1; Mu2 += 2)
				{
					for (int Ml2 = - Jl2_; Ml2 < Jl2_ + 1; Ml2 += 2)
					{										
						if (std::abs(Mu2 - Ml2) <= 2) 
						{
							const int q2 = Ml2 - Mu2;

				      		const double W3J1 = W3JS(Ju2_, Jl2_, 2,-Mu2, Ml2, -q2);  
				      		const double W3J2 = W3JS(2, 2, 2 * K, q2, -q2, 0); 

							const double fact = coeff_K * std::pow(-1.0, double(q2) / 2.0 + 1.0) * std::pow(W3J1, 2) * W3J2;								   
							
		      				const double um = coeff2 * (gu_ * (double(Mu2) / 2.0) - gl_ * (double(Ml2) / 2.0)) + u_red;		      											      					      			

							for (int Q = -K; Q < K + 1; ++Q)
							{			
								const std::complex<double> faddeva = Faddeeva::w(um + a_damp);
				        		const auto D_KQQ                   = std::conj(R(K, 0, Q));

				        		const auto fact_re = fact * std::real(faddeva) * D_KQQ;
				        		const auto fact_im = fact * std::imag(faddeva) * D_KQQ;	

				        		for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
								{

									auto TKQi = get_TKQi(T_KQ[i_stokes], K, Q);										

									// etas
									block_eta[b + i_stokes] += std::real(fact_re * TKQi);										

									// rhos
									if (i_stokes > 0) block_rho[b + i_stokes] += std::real(fact_im * TKQi);
								}						        								
							}
						}
					}		
				}
        	}

        	if (enable_continuum_) block_eta[b] += k_c[n_nu];      

        	// if (i == 0 and j == 0) 
        	// {
        	// 	// std::cout << "k = "     <<   k      << std::endl; 
        	// 	std::cout << block_eta[b] << std::endl;         		
        	// }	         	
        	
        	// sanity checks
        	if (block_eta[b] == 0) std::cerr << "\nWARNING: zero eta_I!"     << std::endl; 
        	if (block_eta[b] < 0)  std::cerr << "\nWARNING: negative eta_I!" << std::endl; 		
        	if (block_rho[b] < 0)  std::cerr << "\nWARNING: negative rho_I!" << std::endl; 	      

        	if (isnan(block_eta[b    ])) std::cerr << "\nWARNING: eta_I = NaN!" << std::endl; 
        	if (isnan(block_eta[b + 1])) std::cerr << "\nWARNING: eta_Q = NaN!" << std::endl; 
        	if (isnan(block_eta[b + 2])) std::cerr << "\nWARNING: eta_U = NaN!" << std::endl; 
        	if (isnan(block_eta[b + 3])) std::cerr << "\nWARNING: eta_V = NaN!" << std::endl;         	
        	if (isnan(block_rho[b + 1])) std::cerr << "\nWARNING: rho_Q = NaN!" << std::endl; 
        	if (isnan(block_rho[b + 2])) std::cerr << "\nWARNING: rho_U = NaN!" << std::endl; 
        	if (isnan(block_rho[b + 3])) std::cerr << "\nWARNING: rho_V = NaN!" << std::endl;         	      	        	
        } 	
    });	

	// debug			
	// const Real dichroism_module = std::sqrt(etas_and_rhos[1] * etas_and_rhos[1] + etas_and_rhos[2] * etas_and_rhos[2] + etas_and_rhos[3] * etas_and_rhos[3]);	
	// if (etas_and_rhos[0] < dichroism_module) dichroism_warning = true;				
	// if (dichroism_warning) std::cout << "\nWARNING: eta_I < eta! (Eq. (7) Gioele Paganini 2018, Part III)" << std::endl; 		
}


void RT_problem::set_TKQ_tensor()
{
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < N_theta_; ++j)
		{
			for (int k = 0; k < N_chi_; ++k)
			{
				auto T_KQ = compute_T_KQ(i, theta_grid_[j], chi_grid_[k]);
			
				T_KQ_.push_back(T_KQ);							
			}
		}
	}
}


void RT_problem::set_up(){
 
    if (mpi_rank_ == 0) std::cout << "\nPrecomputing quantities...";				

    // temporary constants
    Real tmp_const, tmp_const2, tmp_const3;
    
	// compute line-center frequency
	const Real dE = Eu_ - El_;
	nu_0_ = dE * c_;

	// compute coefficients depending on the spatial point 
	
	// const for k_L
	tmp_const = c_ * c_* (Ju2_ + 1) * Aul_ / (8 * PI * nu_0_ * nu_0_ * (Jl2_ + 1));

	// const for Wien 
	tmp_const2 = 2 * h_ * nu_0_ * nu_0_ * nu_0_ / (c_ * c_);
	
	// const for D1
	tmp_const3 = (Ju2_ * Ju2_ + 2 * Ju2_ - 3) / ( 3 * (Ju2_ * Ju2_ + 2 * Ju2_ - 7));

	const Real mass_real = mass_ * 1.6605e-24;

	auto xi_dev   = xi_  ->view_device();
	auto T_dev    = T_   ->view_device();
	auto Nl_dev   = Nl_  ->view_device();
	auto Cul_dev  = Cul_ ->view_device();
	auto Qel_dev  = Qel_ ->view_device();
	auto D2_dev   = D2_  ->view_device();
	auto D1_dev   = D1_  ->view_device();
	auto k_L_dev  = k_L_ ->view_device();
	auto u_dev    = u_   ->view_device();
	auto a_dev    = a_   ->view_device();

	auto Doppler_width_dev = Doppler_width_->view_device();
	auto epsilon_dev       = epsilon_ ->view_device();
	auto W_T_dev           = W_T_->view_device();

	// compute atmospheric quantities 
    sgrid::parallel_for("INIT-ATM", space_grid_->md_range(), KOKKOS_LAMBDA(int i, int j, int k) 
    {       	
    	auto *u = u_dev.block(i, j, k);

    	// assign some variables for readability
    	Real T   =   T_dev.ref(i,j,k);    	    	
    	Real xi  =  xi_dev.ref(i,j,k);
    	Real Cul = Cul_dev.ref(i,j,k);
    	
        // precompute quantities depening only on position
        if (not use_PORTA_input_) 
        {
        	D2_dev.ref(i,j,k)  = 0.5 * Qel_dev.ref(i,j,k); 
        	epsilon_dev.ref(i,j,k) = Cul/(Cul + Aul_);	
        }

		D1_dev.ref(i,j,k)  = tmp_const3 * D2_dev.ref(i,j,k);

		k_L_dev.ref(i,j,k) = tmp_const * Nl_dev.ref(i,j,k);		
		
		Doppler_width_dev.ref(i,j,k) = dE * std::sqrt(xi * xi + 2 * k_B_ * T / mass_real);			

		if (use_PORTA_input_) a_dev.ref(i,j,k) = (Aul_ + Cul + Qel_dev.ref(i,j,k)) / (4 * PI * Doppler_width_dev.ref(i,j,k));

		W_T_dev.ref(i,j,k) = tmp_const2 * std::exp(- h_ * nu_0_ / (k_B_ * T));		
		
		// on position and frequency
		for (int n = 0; n < N_nu_; ++n)
		{
			u[n] = (nu_0_ - nu_grid_[n]) / Doppler_width_dev.ref(i,j,k);						
		}			
    });			      

	// precompute polarization tensors T_KQ
	set_TKQ_tensor();
	   	
	// compute etas and rhos
	set_eta_and_rhos();

	// delete stuff that is not needed -----------------------------------> TODO?
	
	if (mpi_rank_ == 0) std::cout << "done" << std::endl;	
}


void const RT_problem::print_surface_profile(const Field_ptr_t field, const int i_stoke, const int i_space, const int j_space, const int j_theta, const int k_chi){
		
	MPI_Barrier(MPI_COMM_WORLD);
		
	const auto f_dev = field->view_device();	
	const auto g_dev = space_grid_->view_device();

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
				
	if (g_dev.global_coord(2, k_start) == 0)
	{
		for (int i = i_start; i < i_end; ++i)
		{			
			if (g_dev.global_coord(0, i) == i_space)
			{
				for (int j = j_start; j < j_end; ++j)
				{										
					if (g_dev.global_coord(1, j) == j_space)
					{				
						// print info	
						switch (i_stoke)
						{
							case (0): 
								std::cout << "\nSurface radiation, Stoke parameter I, ";
								break;
							case (1): 
								std::cout << "\nSurface radiation, Stoke parameter Q, ";
								break;
							case (2): 
								std::cout << "\nSurface radiation, Stoke parameter U, ";
								break;
							case (3): 
								std::cout << "\nSurface radiation, Stoke parameter V, ";
								break;
							default:
								std::cout << "\nERROR: i_stoke should be smaller then 4!" << std::endl;
						}

						// print info
						std::cout << "i = " << i_space << ", j = " << j_space << ", mu =  " 
								  << mu_grid_[j_theta] << ", chi =  " << chi_grid_[k_chi] 
								  << ", mpi_rank = " << mpi_rank_ << std::endl;	
						
						const int b_start = i_stoke + local_to_block(j_theta, k_chi, 0);						

						for (int b = 0; b < 4 * N_nu_; b = b + 4) 
						{								
							std::cout << f_dev.block(i,j,k_start)[b_start + b] << std::endl; 							
						}
					}
				}
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);	
	std::this_thread::sleep_for(std::chrono::seconds(1));
}


void const RT_problem::print_surface_QI_profile(const Field_ptr_t field, const int i_space, const int j_space, 
	 const int j_theta, const int k_chi, const int i_stokes, const bool center_line){

	MPI_Barrier(MPI_COMM_WORLD);

	if (i_stokes > 3 or i_stokes < 1) std::cout << "ERROR in print_surface_QI_profile input!" << std::endl; 

	if (mpi_rank_ == 0 and i_stokes == 1) std::cout << "\nSurface Q/I, ";
	if (mpi_rank_ == 0 and i_stokes == 2) std::cout << "\nSurface U/I, ";
	if (mpi_rank_ == 0 and i_stokes == 3) std::cout << "\nSurface V/I, ";

	if (mpi_rank_ == 0) std::cout << "mu =  " << mu_grid_[j_theta] << ", chi =  " << chi_grid_[k_chi] << std::endl;		
	if (mpi_rank_ == 0) std::cout << "i,j =  " << i_space << ", " << j_space << std::endl;	
	
	const auto f_dev = field->view_device();	
	const auto g_dev = space_grid_->view_device();

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
		
	int i_global, j_global;
	
	if (g_dev.global_coord(2, k_start) == 0)
	{
		for (int i = i_start; i < i_end; ++i)
		{
			i_global = g_dev.global_coord(0, i);

			if (i_global == i_space)
			{
				for (int j = j_start; j < j_end; ++j)
				{
					j_global = g_dev.global_coord(1, j);

					if (j_global == j_space)
					{
						const int b_start = local_to_block(j_theta, k_chi, 0);						

						if (center_line)
						{
							// print only center line
							const int b_nu_0 = 4 * (N_nu_ - 1) / 2;
							const double I   = f_dev.block(i,j,k_start)[b_start + b_nu_0];
							const double QUV = f_dev.block(i,j,k_start)[b_start + b_nu_0 + i_stokes];	

							std::cout << QUV/I << std::endl; 						
						}
						else
						{
							for (int b = 0; b < 4 * N_nu_; b = b + 4) 
							{															
								const double I   = f_dev.block(i,j,k_start)[b_start + b];
								const double QUV = f_dev.block(i,j,k_start)[b_start + b + i_stokes];							

								std::cout << QUV/I << std::endl; 														
							}
						}						

						std::this_thread::sleep_for(std::chrono::seconds(1));
						MPI_Barrier(MPI_COMM_WORLD);						

						return;
					}
				}
			}
		}
	}

	std::this_thread::sleep_for(std::chrono::seconds(1));
	MPI_Barrier(MPI_COMM_WORLD);	
}



void const RT_problem::print_surface_QI_point(const int i_space, const int j_space, const int j_theta, 
											  const int k_chi, const int n_nu, const int i_stokes)
{
	MPI_Barrier(MPI_COMM_WORLD);

	if (i_stokes > 3 or i_stokes < 1) std::cout << "ERROR in print_surface_QI_point input!" << std::endl; 

	if (mpi_rank_ == 0 and i_stokes == 1) std::cout << "\nSurface Q/I, ";
	if (mpi_rank_ == 0 and i_stokes == 2) std::cout << "\nSurface U/I, ";
	if (mpi_rank_ == 0 and i_stokes == 3) std::cout << "\nSurface V/I, ";

	if (mpi_rank_ == 0) std::cout << "mu =  " << mu_grid_[j_theta] << ", chi =  " << chi_grid_[k_chi] << ", nu =  " << nu_grid_[n_nu] << std::endl;		
	
	const auto f_dev = I_field_->view_device();	
	const auto g_dev = space_grid_->view_device();

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
		
	int i_global, j_global;
	
	if (g_dev.global_coord(2, k_start) == 0)
	{
		for (int i = i_start; i < i_end; ++i)
		{
			i_global = g_dev.global_coord(0, i);

			if (i_global == i_space)
			{
				for (int j = j_start; j < j_end; ++j)
				{
					j_global = g_dev.global_coord(1, j);

					if (j_global == j_space)
					{
						const int b_start = local_to_block(j_theta, k_chi, n_nu);
						const double I   = f_dev.block(i,j,k_start)[b_start];
						const double QUV = f_dev.block(i,j,k_start)[b_start + i_stokes];							

						std::cout << QUV/I << std::endl; 													

						return;
					}
				}
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	std::this_thread::sleep_for(std::chrono::seconds(1));
}



void const RT_problem::print_profile(const Field_ptr_t field, const int i_stoke, 
									 const int i_space, const int j_space, const int k_space, 
									 const int j_theta, const int k_chi)
{

	if (mpi_rank_ == 0)
	{
		switch (i_stoke)
		{
			case (0): 
				std::cout << "\nSurface radiation, Stoke parameter I, ";
				break;
			case (1): 
				std::cout << "\nSurface radiation, Stoke parameter Q, ";
				break;
			case (2): 
				std::cout << "\nSurface radiation, Stoke parameter U, ";
				break;
			case (3): 
				std::cout << "\nSurface radiation, Stoke parameter V, ";
				break;
			default:
				std::cout << "\nERROR: i_stoke should be smaller then 4!" << std::endl;
		}

		std::cout << "mu =  "       << mu_grid_[j_theta];
		std::cout << ", chi =  "    << chi_grid_[k_chi];		
		std::cout << ", height =  " << depth_grid_[k_space] << " km" << std::endl;				
	}

	const auto f_dev = field->view_device();	
	const auto g_dev = space_grid_->view_device();

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
	const int k_end = k_start + g_dev.dim[2];
		
	int i_global, j_global, k_global;

	for (int i = i_start; i < i_end; ++i)
	{
		i_global = g_dev.global_coord(0, i);

		if (i_global == i_space)
		{
			for (int j = j_start; j < j_end; ++j)
			{
				j_global = g_dev.global_coord(1, j);

				if (j_global == j_space)
				{
					for (int k = k_start; k < k_end; ++k)
					{
						k_global = g_dev.global_coord(2, k);

						if (k_global == k_space)
						{
							const int b_start = i_stoke + local_to_block(j_theta, k_chi, 0);

							for (int b = 0; b < 4 * N_nu_; b = b + 4) 
							{	
								std::cout << f_dev.block(i,j,k)[b_start + b] << std::endl; 							
							}					
						}
					}					

					return;
				}
			}
		}
	}	

	MPI_Barrier(MPI_COMM_WORLD);
}

// write surface profile in one single point
void const RT_problem::write_surface_point_profiles(input_string file_name, const int i_space, const int j_space)
{
	// // a single MPI rank writes output
	// if (mpi_rank_ == 0) std::cout << " Writing output in spatial point (" << i_space << ", " << j_space << ")" << std::endl;

	const auto f_dev = I_field_->view_device();	
	const auto g_dev = space_grid_->view_device();

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
		
	int i_global, j_global;

	double I, QUV;

	// write profiles
	if (g_dev.global_coord(2, k_start) == 0)
	{
		for (int i = i_start; i < i_end; ++i)
		{
			i_global = g_dev.global_coord(0, i);

			if (i_global == i_space)
			{
				for (int j = j_start; j < j_end; ++j)
				{
					j_global = g_dev.global_coord(1, j);

					if (j_global == j_space)
					{						
						// Create a new file 
						input_string output_file = file_name + "_" + std::to_string(i_space) + "_" + std::to_string(j_space) + ".m";
					 	std::ofstream outputFile(output_file);

						if (outputFile.is_open()) 		
						{		
							// write grids 
							outputFile << "\nnu_grid_ = [ ";
							for (int n = 0; n < N_nu_; ++n) outputFile << std::scientific << std::setprecision(15) << nu_grid_[n] << " ";
							outputFile << "];\n";

							outputFile << "\ntheta_grid = [ ";
							for (int j_theta = 0; j_theta < N_theta_; ++j_theta) outputFile << std::scientific << std::setprecision(15) << theta_grid_[j_theta] << " ";
							outputFile << "];\n";

							outputFile << "\nmu_grid = [ ";
							for (int j_theta = 0; j_theta < N_theta_; ++j_theta) outputFile << std::scientific << std::setprecision(15) << mu_grid_[j_theta] << " ";
							outputFile << "];\n";

							outputFile << "\nchi_grid = [ ";
							for (int k_chi = 0; k_chi < N_chi_; ++k_chi) outputFile << std::scientific << std::setprecision(15) << chi_grid_[k_chi] << " ";
							outputFile << "];\n";

							// create MATLAB data structure
							outputFile << std::scientific << std::setprecision(15) <<  "\nField = cell(4," << N_theta_ << "," << N_chi_ << ");" << std::endl;

							for (int j_theta = N_theta_/2; j_theta < N_theta_; ++j_theta)
							{								
								for (int k_chi = 0; k_chi < N_chi_; ++k_chi)
								{							
									const int b_start = local_to_block(j_theta, k_chi, 0);						
						
									for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
									{
										outputFile << "\nField{" << i_stokes + 1 << "," << j_theta + 1 << "," << k_chi + 1 << "} = [ ";										

										for (int b = 0; b < 4 * N_nu_; b = b + 4) 
										{
											I = f_dev.block(i,j,k_start)[b_start + b];

											if (i_stokes == 0)
											{
												outputFile << std::scientific << std::setprecision(15) << I << " ";
											}
											else
											{
												QUV = f_dev.block(i,j,k_start)[b_start + b + i_stokes];

												outputFile << std::scientific << std::setprecision(15) << 100.0 * QUV/I << " ";
											}
										}

										outputFile << "];\n";	
									}
								}
							}

							// Close the file		
    						outputFile.close();

    						if (mpi_rank_ == 0) std::cout << "Output written in " << output_file << "\n" << std::endl;    						
    					} 
						else
					  	{
					    	if (mpi_rank_ == 0) std::cout << "\nERROR: failed to create the output file." << std::endl;
					  	}
					}
				}
			}
		}
	}  	
}


// write surface profile in all surface
void const RT_problem::write_surface_profiles(input_string file_name)
{	
	if (mpi_size_x_ > 1 or mpi_size_y_ > 1)
	{
		std::cerr << "\nWARNING: write_surface_profiles not supported for hotizontal decomposition!" << std::endl;    						
	}
	
	const auto f_dev = I_field_->view_device();	
	const auto g_dev = space_grid_->view_device();

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
		
	int i_global, j_global;

	double I, QUV;

	// write profiles
	if (g_dev.global_coord(2, k_start) == 0)
	{
		// Create a new file 
		input_string output_file = file_name + ".m";
		std::ofstream outputFile(output_file);

		if (outputFile.is_open()) 		
		{		
			if (mpi_rank_ == 0) std::cout << "\nWriting surface output in " << output_file << "\n" << std::endl;    						

			// write grids 
			outputFile << "\nnu_grid_ = [ ";
			for (int n = 0; n < N_nu_; ++n) outputFile << nu_grid_[n] << " ";
			outputFile << "];\n";

			outputFile << "\ntheta_grid = [ ";
			for (int j_theta = 0; j_theta < N_theta_; ++j_theta) outputFile << theta_grid_[j_theta] << " ";
			outputFile << "];\n";

			outputFile << "\nmu_grid = [ ";
			for (int j_theta = 0; j_theta < N_theta_; ++j_theta) outputFile << mu_grid_[j_theta] << " ";
			outputFile << "];\n";

			outputFile << "\nchi_grid = [ ";
			for (int k_chi = 0; k_chi < N_chi_; ++k_chi) outputFile << chi_grid_[k_chi] << " ";
			outputFile << "];\n";

			// create MATLAB data structure
			outputFile <<  "\nField = cell(4," << N_x_ << "," << N_y_ << "," << N_theta_ << "," << N_chi_ << ");" << std::endl;			
		}
		else
		{
			if (mpi_rank_ == 0) std::cout << "\nERROR: failed to create the output file." << std::endl;
		}

		// loop over spatial ppoints and directions
		for (int i = i_start; i < i_end; ++i)
		{
			i_global = g_dev.global_coord(0, i);
			
			for (int j = j_start; j < j_end; ++j)
			{
				j_global = g_dev.global_coord(1, j);
											
				for (int j_theta = N_theta_/2; j_theta < N_theta_; ++j_theta)
				{								
					for (int k_chi = 0; k_chi < N_chi_; ++k_chi)
					{							
						const int b_start = local_to_block(j_theta, k_chi, 0);						
				
						for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
						{
							outputFile << "\nField{" << i_stokes + 1 << "," << i_global + 1 << "," << j_global + 1<< "," << j_theta + 1 << "," << k_chi + 1 << "} = [ ";

							for (int b = 0; b < 4 * N_nu_; b = b + 4) 
							{
								I = f_dev.block(i,j,k_start)[b_start + b];

								if (i_stokes == 0)
								{
									outputFile << I << " ";
								}
								else
								{
									QUV = f_dev.block(i,j,k_start)[b_start + b + i_stokes];

									outputFile << 100.0 * QUV/I << " ";
								}
							}

							outputFile << "];\n";	
						}
					}
				}				
			}		
		}
		
		// Close the file		
		outputFile.close();
	}  		
}


// write surface profile in one single point
void const RT_problem::write_surface_point_profiles_Omega(input_string file_name, const int i_space, const int j_space)
{
	// // a single MPI rank writes output
	// if (mpi_rank_ == 0) std::cout << " Writing output in spatial point (" << i_space << ", " << j_space << ")" << std::endl;

	const auto I_dev = I_field_Omega_->view_device();	
	const auto g_dev = space_grid_->view_device();

	const int block_size = 4 * N_nu_;

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
		
	int i_global, j_global;

	double I, QUV;

	// write profiles
	if (g_dev.global_coord(2, k_start) == 0)
	{
		for (int i = i_start; i < i_end; ++i)
		{
			i_global = g_dev.global_coord(0, i);

			if (i_global == i_space)
			{
				for (int j = j_start; j < j_end; ++j)
				{
					j_global = g_dev.global_coord(1, j);

					if (j_global == j_space)
					{						
						// Create a new file 
						input_string output_file = file_name + "_" + std::to_string(i_space) + "_" + std::to_string(j_space) + ".m";
					 	std::ofstream outputFile(output_file);

						if (outputFile.is_open()) 		
						{									
							// create MATLAB data structure							
							outputFile <<  "\nField_Omega = cell(4,1);" << std::endl;
						
								for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
								{
									outputFile << "\nField_Omega{" << i_stokes + 1 << "} = [ ";										

									for (int b = 0; b < block_size; b = b + 4) 
									{
										I = I_dev.block(i,j,k_start)[b];

										if (i_stokes == 0)
										{
											outputFile << I << " ";
										}
										else
										{
											QUV = I_dev.block(i,j,k_start)[b + i_stokes];

											outputFile << 100.0 * QUV/I << " ";
										}
									}

									outputFile << "];\n";	
								}													

							// Close the file		
    						outputFile.close();

    						if (mpi_rank_ == 0) std::cout << "Output written in " << output_file << "\n" << std::endl;    						
    					} 
						else
					  	{
					    	if (mpi_rank_ == 0) std::cout << "\nERROR: failed to create the output file." << std::endl;
					  	}
					}
				}
			}
		}
	}  	
}


// write surface profile for arbitrary direction
void const RT_problem::write_surface_profiles_Omega(input_string file_name)
{	
	if (mpi_size_x_ > 1 or mpi_size_y_ > 1)
	{
		std::cerr << "\nWARNING: write_surface_profiles_Omega not supported for hotizontal decomposition!" << std::endl;    						
	}

	const auto I_dev = I_field_Omega_->view_device();	
	const auto g_dev = space_grid_->view_device();

	const int block_size = 4 * N_nu_;

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
		
	int i_global, j_global;

	double I, QUV;

	// write profiles
	if (g_dev.global_coord(2, k_start) == 0)
	{
		// Create a new file 
		input_string output_file = file_name + ".m";
		std::ofstream outputFile(output_file);

		if (outputFile.is_open()) 							
		{						
			if (mpi_rank_ == 0) std::cout << "\nWriting surface output in " << output_file << "\n" << std::endl;    						
			
			// create MATLAB data structure
			outputFile <<  "\nField = cell(4," << N_x_ << "," << N_y_ << ");" << std::endl;					
		}
		else
		{
			if (mpi_rank_ == 0) std::cout << "\nERROR: failed to create the output file." << std::endl;
		}	

		for (int i = i_start; i < i_end; ++i)
		{
			i_global = g_dev.global_coord(0, i);

			for (int j = j_start; j < j_end; ++j)
			{
				j_global = g_dev.global_coord(1, j);
																							
				for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
				{
					outputFile << "\nField{" << i_stokes + 1 << "," << i_global + 1 << "," << j_global + 1 << "} = [ ";

					for (int b = 0; b < block_size; b = b + 4) 
					{
						I = I_dev.block(i,j,k_start)[b];

						if (i_stokes == 0)
						{
							outputFile << I << " ";
						}
						else
						{
							QUV = I_dev.block(i,j,k_start)[b + i_stokes];

							outputFile << 100.0 * QUV/I << " ";
						}
					}

					outputFile << "];\n";	
				}																									
			}			
		}

		// Close the file		
    	outputFile.close();
	}  	
}




bool RT_problem::field_is_zero(const Field_ptr_t field)
{
	bool field_is_zero = true;

	auto g_dev = space_grid_->view_device();
	auto f_dev =       field->view_device();	

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
	const int k_end = k_start + g_dev.dim[2];	
	
	for (int k = k_start; k < k_end; ++k)					
	{															
		for (int j = j_start; j < j_end; ++j)
		{
			for (int i = i_start; i < i_end; ++i)				
			{
				for (int b = 0; b < block_size_; b++) 
				{
					if (f_dev.block(i,j,k)[b] != 0)
					{
						field_is_zero = false;
						break;
					}
				}							
			}
		}
	}

    return field_is_zero;
}


void RT_problem::set_grid_partition() // TODO remove hardcoding
{
	// TODO: now hardcoded	

	if (mpi_size_ <= N_z_)
	{		
		mpi_size_z_ = mpi_size_;
		mpi_size_x_ = 1;
		mpi_size_y_ = 1;
	}
	else
	{			
		if (mpi_rank_ == 0) std::cout << "========== WARNING: hardcoding grid partition ==========" << std::endl;					

		if (mpi_size_ == 1152) // HARDCODED
		{
			mpi_size_z_ = 128;
			mpi_size_x_ = 3;
			mpi_size_y_ = 3;			
		}	
		else if (mpi_size_ == 384) // HARDCODED
		{
			mpi_size_z_ = 128;
			mpi_size_x_ = 3;
			mpi_size_y_ = 1;			
		}
		else if (mpi_size_ == 256) // HARDCODED
		{
			mpi_size_z_ = 128;
			mpi_size_x_ = 2;
			mpi_size_y_ = 1;			
		}
		else if (mpi_size_ == 512) // HARDCODED
		{
			mpi_size_z_ = 128;
			mpi_size_x_ = 2;
			mpi_size_y_ = 2;			
		}
		else if (mpi_size_ == 1024) // HARDCODED
		{
			mpi_size_z_ = 128;
			mpi_size_x_ = 4;
			mpi_size_y_ = 2;	
		}
		else if (mpi_size_ == 2048) // HARDCODED
		{
			mpi_size_z_ = 128;
			mpi_size_x_ = 4;
			mpi_size_y_ = 4;		
		}
		else if (mpi_size_ == 4096) // HARDCODED
		{
			mpi_size_z_ = 1;
			mpi_size_x_ = 64;
			mpi_size_y_ = 64;			
		}
		else if (mpi_size_ == 6144) // HARDCODED
		{
			mpi_size_z_ = 128;
			mpi_size_x_ = 4;
			mpi_size_y_ = 12;			
		}
		else if (mpi_size_ == 8192) // HARDCODED
		{
			mpi_size_z_ = 128;
			mpi_size_x_ = 8;
			mpi_size_y_ = 8;			
		}
		else if (mpi_size_ == 12288) // HARDCODED
		{
			mpi_size_z_ = 128;
			mpi_size_x_ = 12;
			mpi_size_y_ = 8;			
		}
		else if (mpi_size_ == 16384) 
		{
			mpi_size_z_ = 128;
            mpi_size_x_ = 16;
			mpi_size_y_ = 8;
		}
		else
		{
			mpi_size_z_ = N_z_;

			const double mpi_size_xy = mpi_size_ / mpi_size_z_;

			if (std::floor(mpi_size_xy) != mpi_size_xy) std::cout << "ERROR: problem with domain decomposition: mpi_size_ / mpi_size_z_ not integer!" << std::endl;

			mpi_size_x_ = std::floor(std::sqrt(mpi_size_xy));

			const double mpi_size_y = mpi_size_xy / mpi_size_x_;

			if (std::floor(mpi_size_y) != mpi_size_y) std::cout << "ERROR: problem with domain decomposition: mpi_size_y not integer!" << std::endl;

			mpi_size_y_ = (int) mpi_size_y;		
		}
	}
}




