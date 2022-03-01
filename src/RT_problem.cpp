#include "RT_problem.hpp"


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

	Jl_ = Jl2_/2;
	Ju_ = Ju2_/2;
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

	const size_t N_z_N_nu_ = N_z_ * N_nu_;

	// safety check
	if (sigma_vec.size()    != N_z_N_nu_) std::cout << "WARNING: size mismatch in read_continumm_1D()" << std::endl;
	if (k_c_vec.size()      != N_z_N_nu_) std::cout << "WARNING: size mismatch in read_continumm_1D()" << std::endl;
	if (eps_c_th_vec.size() != N_z_N_nu_) std::cout << "WARNING: size mismatch in read_continumm_1D()" << std::endl;
	
	auto g_dev = space_grid_->view_device();

	// fill field
	sgrid::parallel_for("READ SIGMA", space_grid_->md_range(), SGRID_LAMBDA(int i, int j, int k) {

		int k_global;

		for (int n = 0; n < (int)N_nu_; ++n)
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
			entry *= 1399600; // convert to Larmor frequency
			nu_L_vec.push_back(entry);
			lineStream >> entry;
			theta_B_vec.push_back(entry);
			lineStream >> entry;
			chi_B_vec.push_back(entry);				
		}		

		first_line = false;
	} 

	// safety check
	if (nu_L_vec.size()    != N_z_) std::cout << "WARNING: size mismatch in read_magnetic_field_1D()" << std::endl;
	if (theta_B_vec.size() != N_z_) std::cout << "WARNING: size mismatch in read_magnetic_field_1D()" << std::endl;
	if (chi_B_vec.size()   != N_z_) std::cout << "WARNING: size mismatch in read_magnetic_field_1D()" << std::endl;
	
	auto g_dev = space_grid_->view_device();

	// fill field
	sgrid::parallel_for("READ B", space_grid_->md_range(), SGRID_LAMBDA(int i, int j, int k) {
				
		const int k_global = g_dev.start[2] + k - g_dev.margin[2];

		B_dev.block(i, j, k)[0] =    nu_L_vec[k_global];					
		B_dev.block(i, j, k)[1] = theta_B_vec[k_global];					
		B_dev.block(i, j, k)[2] =   chi_B_vec[k_global];													
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

	// fill field 
	sgrid::parallel_for("READ-ATM1D", space_grid_->md_range(), SGRID_LAMBDA(int i, int j, int k) {

		const int k_global = g_dev.start[2] + k - g_dev.margin[2];
								
		T_dev.ref(  i, j, k) =   T_vec[k_global];
		xi_dev.ref( i, j, k) =  xi_vec[k_global];		
		a_dev.ref(  i, j, k) =   a_vec[k_global];		
		Nl_dev.ref( i, j, k) =  Nl_vec[k_global];		
		Cul_dev.ref(i, j, k) = Cul_vec[k_global];		
		Qel_dev.ref(i, j, k) = Qel_vec[k_global];						
	});
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


void RT_problem::read_frequency(input_string filename){

	if (mpi_rank_ == 0) std::cout << "Reading frequencies [s-1] from " << filename << std::endl;

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
			lineStream >> entry;			
			lineStream >> entry;
			nu_grid_.push_back(entry);		

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

	block_size_ = 4 * N_nu_ * N_theta_ * N_chi_;
	tot_size_   = N_s_ * block_size_;	

	if (mpi_rank_ == 0 and mpi_size_ > (int) N_s_) std::cerr << "\n========= WARNING: mpi_size > N_s! =========\n" << std::endl;
}

void const RT_problem::print_info(){
	
	if (mpi_rank_ == 0) 		
	{		
		std::cout << "\n=========== 2-levels atom parameters ===========\n" << std::endl;		
		std::cout << "Mass = " << mass_ << std::endl;
		std::cout << "El = "   << El_   << std::endl;
		std::cout << "Eu = "   << Eu_   << std::endl;		
		std::cout << "Jl = "   << Jl_   << std::endl;
		std::cout << "Ju = "   << Ju_   << std::endl;
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
		
		std::cout << "\ntotal size = " << tot_size_   << std::endl;			
		std::cout << "block size = "   << block_size_ << std::endl;		

		std::cout << "\ntheta grid = [ ";

		for (int i = 0; i < (int)N_theta_; ++i) std::cout << theta_grid_[i] << " ";

		std::cout << "]\nmu grid    = [ ";

		for (int i = 0; i < (int)N_theta_; ++i) std::cout   << mu_grid_[i] << " ";

		std::cout << "]\nchi grid   = [ ";

		for (int i = 0; i < (int)N_chi_; ++i) std::cout   << chi_grid_[i] << " ";

		std::cout << "] " << std::endl;		

		if (only_vertical_decomposition_) std::cout << "\nDomain decompostion only in the z direction (Jiri method)" << std::endl;		
	}
}


void RT_problem::allocate_fields(){

		// create fields 
		I_field_ = std::make_shared<Field_t>("I", space_grid_, block_size_, sgrid::BOX_STENCIL);
		S_field_ = std::make_shared<Field_t>("S", space_grid_, block_size_, sgrid::BOX_STENCIL);

		eta_field_ = std::make_shared<Field_t>("eta", space_grid_, block_size_, sgrid::BOX_STENCIL);
		rho_field_ = std::make_shared<Field_t>("rho", space_grid_, block_size_, sgrid::BOX_STENCIL);

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
	B_ = std::make_shared<Field_t>("B", space_grid_, 3, sgrid::STAR_STENCIL); 	

	// bulk velocities, in polar coordinates
	v_b_ = std::make_shared<Field_t>("v_b", space_grid_, 3, sgrid::STAR_STENCIL);    
	
	// quantities depending on position that can be precomputed
	Doppler_width_ = std::make_shared<Field_t>("Doppler_width", space_grid_);
	k_L_           = std::make_shared<Field_t>("k_L", 		    space_grid_);
	epsilon_       = std::make_shared<Field_t>("epsilon_", 		space_grid_);

	// input quantities depending on position and frequency 
	u_        = std::make_shared<Field_t>("u",        space_grid_, N_nu_, sgrid::STAR_STENCIL);  
	sigma_    = std::make_shared<Field_t>("sigma",    space_grid_, N_nu_, sgrid::STAR_STENCIL);
	k_c_      = std::make_shared<Field_t>("k_c",      space_grid_, N_nu_, sgrid::STAR_STENCIL);
	eps_c_th_ = std::make_shared<Field_t>("eps_c_th", space_grid_, N_nu_, sgrid::STAR_STENCIL);
	
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

    Kokkos::parallel_for("INIT I", space_grid_->md_range(), KOKKOS_LAMBDA(int i, int j, int k) 
    {         
        auto *block = field_dev.block(i, j, k);
         
        for (int b = 0; b < (int)block_size_; ++b) 
        {
        	block[b] = input_value;        	
        }
    });
}


void RT_problem::set_theta_chi_grids(const size_t N_theta, const size_t N_chi, const bool double_GL){

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
    	for (size_t i = 0; i < N_theta / 2; ++i)
	    {
	    	mu_grid_.push_back(mu_1[i]);
	    	w_theta_.push_back(w_1[i]);

			theta_grid_.push_back(std::acos(mu_1[i]));		    	
	    }

	    // second half
    	for (size_t i = 0; i < N_theta / 2; ++i)
	    {
	    	mu_grid_.push_back(mu_2[i]);
	    	w_theta_.push_back(w_2[i]);

			theta_grid_.push_back(std::acos(mu_2[i]));		    	
	    }
    }
    else
    {
    	legendre_rule(N_theta,  -1.0, 1.0, mu_grid_, w_theta_);

	    for (size_t i = 0; i < N_theta; ++i)
	    {
	    	theta_grid_.push_back(std::acos(mu_grid_[i]));
	    }
    }

    // init equidistant chi grid in [0, 2pi] and trap weights

    if (N_chi % 2 != 0) std::cout << "WARNING: chi grid is not even!" << std::endl;
    if (N_chi % 4 != 0) std::cout << "WARNING: chi grid is not the same for each quadrant!" << std::endl;
    
    const Real delta_chi = 2.0 * PI / N_chi;

    for (size_t i = 0; i < N_chi; ++i)
    {
    	chi_grid_.push_back((i + 0.5) * delta_chi);	 // avoiding grid axes   	
    	w_chi_.push_back(delta_chi);
    }    
}


std::vector<std::complex<Real> > RT_problem::compute_T_KQ(const size_t stokes_i, const Real theta, const Real chi){

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


std::complex<Real> RT_problem::get_TKQi(const size_t i_stokes, const int K, const int Q, const size_t j, const size_t k){

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

    Kokkos::parallel_for("INIT ETA-RHO", space_grid_->md_range(), KOKKOS_LAMBDA(int i, int j, int k) 
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
        std::vector<size_t> local_idx;
        size_t j_theta, k_chi, n_nu;

		// coeffs
		int q;
		Real fact, coeff_K, W3J1, W3J2, um, v_dot_Omega;
		Real u_red, u_b;

		std::complex<Real> D_KQQ, faddeva, fact_re, fact_im;
         
        for (size_t b = 0; b < block_size_; b = b + 4) 
        {        	
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
			v_dot_Omega = v_b[0] * ( cos(theta_v_b) * cos(theta) + sin(theta_v_b) * sin(theta) * cos(chi - chi_v_b));
			u_b = nu_0_ * v_dot_Omega / (c_ * Doppler_width);

			u_red = u[n_nu] + u_b;			

			for (int K = 0; K < 3; ++K)
			{
				coeff_K = coeff * std::sqrt(3 * ( 2 * K + 1));
	
				for (int Mu2 = - Ju2_; Mu2 < Ju2_ + 1; Mu2 += 2)
				{
					for (int Ml2 = - Jl2_; Ml2 < Jl2_ + 1; Ml2 += 2)
					{				
						const int Mu  = Mu2/2;
						const int Ml  = Ml2/2;

						if (std::abs(Mu - Ml) <= 1)
						{
							q = Ml - Mu;

							W3J1 = W3JS(Ju2_, Jl2_, 2, -Mu2, Ml2, -2*q);
							W3J2 = W3JS(2, 2, 2*K, 2*q, -2*q, 0);

							fact = coeff_K * std::pow(-1,q + 1) * std::pow(W3J1,2) * W3J2;

							um = coeff2 * (gu_ * Mu - gl_ * Ml) + u_red;
					
							for (int Q = -K; Q < K + 1; ++Q)
							{			
								faddeva = Faddeeva::w(um + a_damp);								
								D_KQQ   = std::conj(R(K, 0, Q));

								fact_re = fact * std::real(faddeva) * D_KQQ;
								fact_im = fact * std::imag(faddeva) * D_KQQ;

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

        	// checks
        	if (block_eta[b] == 0) std::cerr << "\nWARNING: zero eta_I!"     << std::endl; 
        	if (block_eta[b] < 0)  std::cerr << "\nWARNING: negative eta_I!" << std::endl; 		
        	if (block_rho[b] < 0)  std::cerr << "\nWARNING: negative rho_I!" << std::endl; 	            	
        } 	
    });	

	// exchange ghost layers (for periodic boundary)
	eta_field_->exchange_halos(); 
	rho_field_->exchange_halos(); 	

	// debug			

	// const Real dichroism_module = std::sqrt(etas_and_rhos[1] * etas_and_rhos[1] + etas_and_rhos[2] * etas_and_rhos[2] + etas_and_rhos[3] * etas_and_rhos[3]);
	
	// if (etas_and_rhos[0] < dichroism_module) dichroism_warning = true;
				
	// if (dichroism_warning) std::cout << "\nWARNING: eta_I < eta! (Eq. (7) Gioele Paganini 2018, Part III)" << std::endl; 	
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
	tmp_const = c_ * c_* (2 * Ju_ + 1) * Aul_ / (8 * PI * nu_0_ * nu_0_ * (2 * Jl_ + 1));

	// const for Wien 
	tmp_const2 = 2 * h_ * nu_0_ * nu_0_ * nu_0_ / (c_ * c_);
	
	// const for D1
	tmp_const3 = (4 * Ju_ * Ju_ + 4 * Ju_ - 3) / ( 3 * (4 * Ju_ * Ju_ + 4 * Ju_ - 7));

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

	auto Doppler_width_dev = Doppler_width_->view_device();
	auto epsilon_dev       = epsilon_ ->view_device();
	auto W_T_dev           = W_T_->view_device();

	// compute atmospheric quantities 
    Kokkos::parallel_for("INIT-ATM", space_grid_->md_range(), KOKKOS_LAMBDA(int i, int j, int k) 
    {       	
    	auto *u = u_dev.block(i, j, k);

    	// assign some variables for readability
    	Real T   =   T_dev.ref(i,j,k);    	    	
    	Real xi  =  xi_dev.ref(i,j,k);
    	Real Cul = Cul_dev.ref(i,j,k);
    	
        // precompute quantities depening only on position
        D2_dev.ref(i,j,k)  = 0.5 * Qel_dev.ref(i,j,k); 

		D1_dev.ref(i,j,k)  = tmp_const3 * D2_dev.ref(i,j,k);

		k_L_dev.ref(i,j,k) = tmp_const * Nl_dev.ref(i,j,k);
		
		Doppler_width_dev.ref(i,j,k) = dE * std::sqrt(2 * k_B_ * T / mass_real + xi * xi);

		epsilon_dev.ref(i,j,k) = Cul/(Cul + Aul_);	

		W_T_dev.ref(i,j,k) = tmp_const2 * std::exp(- h_ * nu_0_ / (k_B_ * T));		

		// on position and frequency
		for (size_t n = 0; n < N_nu_; ++n)
		{
			u[n] = (nu_0_ - nu_grid_[n]) / Doppler_width_dev.ref(i,j,k);			
		}	
    });	
		    
	// precompute polarization tensors T_KQ
	for (int i = 0; i < 4; ++i)
	{
		for (size_t j = 0; j < N_theta_; ++j)
		{
			for (size_t k = 0; k < N_chi_; ++k)
			{
				auto T_KQ = compute_T_KQ(i, theta_grid_[j], chi_grid_[k]);
			
				T_KQ_.push_back(T_KQ);							
			}
		}
	}
	   	
	// compute etas and rhos
	set_eta_and_rhos();

	// delete stuff that is not needed -----------------------------------> TODO?
	
	if (mpi_rank_ == 0) std::cout << "done" << std::endl;	
}






