#include "RT_problem.hpp"


void RT_problem::read_frequency(std::string filename){

	if (mpi_rank_ == 0) std::cout << "Reading frequencies [s-1] from " << filename << std::endl;

	std::ifstream myFile(filename);
	std::string line;	

	if (not myFile.good()) std::cerr << "\nERROR: File " << filename << " does not exist!\n" << std::endl;

	bool first_line = true;

	while(getline(myFile, line))
	{
		std::istringstream lineStream(line);
		double entry;
		std::string line_entry;

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
	// N_s_     = depth_grid_.size();
	N_nu_	 = nu_grid_.size();
	N_theta_ = theta_grid_.size();
	N_chi_	 = chi_grid_.size();
	N_dirs_  = N_theta_ * N_chi_; 

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
}

void RT_problem::allocate_atmosphere(){

	// create spatial coordinates
	space_coord = std::make_shared<Field_t>("xyz",   space_grid_, 3, sgrid::BOX_STENCIL);

	// create atmospheric quantities 
	D1_   = std::make_shared<Field_t>("D1",   space_grid_, 1, sgrid::BOX_STENCIL);
	D2_   = std::make_shared<Field_t>("D2",   space_grid_, 1, sgrid::BOX_STENCIL);
	Nl_   = std::make_shared<Field_t>("Nl",   space_grid_, 1, sgrid::BOX_STENCIL);
	Nu_   = std::make_shared<Field_t>("Nu",   space_grid_, 1, sgrid::BOX_STENCIL);
	T_    = std::make_shared<Field_t>("T",    space_grid_, 1, sgrid::BOX_STENCIL);
	xi_   = std::make_shared<Field_t>("xi",   space_grid_, 1, sgrid::BOX_STENCIL);
	nu_L_ = std::make_shared<Field_t>("nu_L", space_grid_, 1, sgrid::BOX_STENCIL);
	Cul_  = std::make_shared<Field_t>("Cul",  space_grid_, 1, sgrid::BOX_STENCIL);
	Qel_  = std::make_shared<Field_t>("Qel",  space_grid_, 1, sgrid::BOX_STENCIL);
	a_    = std::make_shared<Field_t>("a",    space_grid_, 1, sgrid::BOX_STENCIL);
	W_T_  = std::make_shared<Field_t>("W_T",  space_grid_, 1, sgrid::BOX_STENCIL);

	// magnetic field direction, in polar coordinates   
	theta_B_ = std::make_shared<Field_t>("theta_B", space_grid_, 1, sgrid::BOX_STENCIL); 
	chi_B_   = std::make_shared<Field_t>("chi_B",   space_grid_, 1, sgrid::BOX_STENCIL);

	// bulk velocities, in polar coordinates
	v_b_       = std::make_shared<Field_t>("v_b",       space_grid_, 1, sgrid::BOX_STENCIL);    
	theta_v_b_ = std::make_shared<Field_t>("theta_v_b", space_grid_, 1, sgrid::BOX_STENCIL);
	chi_v_b_   = std::make_shared<Field_t>("chi_v_b",   space_grid_, 1, sgrid::BOX_STENCIL);   
	
	// quantities depending on position that can be precomputed
	Doppler_width_ = std::make_shared<Field_t>("Doppler_width", space_grid_, 1, sgrid::BOX_STENCIL);
	k_L_           = std::make_shared<Field_t>("k_L", 		    space_grid_, 1, sgrid::BOX_STENCIL);
	epsilon_       = std::make_shared<Field_t>("epsilon_", 		space_grid_, 1, sgrid::BOX_STENCIL);

	// input quantities depending on position and frequency 
	u_        = std::make_shared<Field_t>("u",        space_grid_, N_nu_, sgrid::BOX_STENCIL);  
	sigma_    = std::make_shared<Field_t>("sigma",    space_grid_, N_nu_, sgrid::BOX_STENCIL);
	k_c_      = std::make_shared<Field_t>("k_c",      space_grid_, N_nu_, sgrid::BOX_STENCIL);
	eps_c_th_ = std::make_shared<Field_t>("eps_c_th", space_grid_, N_nu_, sgrid::BOX_STENCIL);

	// quantities depending on position and direction
	T_KQ_ = std::make_shared<Field_complex_t>("Doppler_width", space_grid_complex_, 4 * 6 * N_dirs_, sgrid::BOX_STENCIL);  

	// allocate
	D1_  -> allocate_on_device(); 
	D2_  -> allocate_on_device(); 
	Nl_  -> allocate_on_device(); 
	Nu_  -> allocate_on_device(); 
	T_   -> allocate_on_device(); 
	xi_  -> allocate_on_device(); 
	nu_L_-> allocate_on_device(); 
	Cul_ -> allocate_on_device(); 
	Qel_ -> allocate_on_device(); 
	a_   -> allocate_on_device(); 
	W_T_ -> allocate_on_device(); 
	
	theta_B_ -> allocate_on_device(); // TODO, put in one? 
	chi_B_   -> allocate_on_device(); 
	v_b_     -> allocate_on_device(); 
	
	theta_v_b_ -> allocate_on_device(); 
	chi_v_b_   -> allocate_on_device(); 
	
	Doppler_width_ -> allocate_on_device(); 
	k_L_           -> allocate_on_device(); 
	epsilon_       -> allocate_on_device(); 

	u_        -> allocate_on_device(); 
	sigma_ 	  -> allocate_on_device(); 
	k_c_      -> allocate_on_device(); 
	eps_c_th_ -> allocate_on_device(); 
	T_KQ_     -> allocate_on_device(); 
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
    // legendre_rule(N_theta,  0.0, pi_, theta_grid_, w_theta_);

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
    const Real delta_chi = 2.0 * pi_ / N_chi;

    for (size_t i = 0; i < N_chi; ++i)
    {
    	chi_grid_.push_back(i * delta_chi);	    	
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
	complex<Real> fa, fb;

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

void RT_problem::set_up(){
 
    if (mpi_rank_ == 0) std::cout << "\nPrecomputing quantities...";				

    // temporary constants
    Real tmp_const, tmp_const2, tmp_const3;
    
	// compute line-center frequency
	const Real dE = Eu_ - El_;
	nu_0_ = dE * c_;

	// compute coefficients depending on the spatial point //////////////////////////
	
	// const for k_L
	tmp_const = c_ * c_* (2 * Ju_ + 1) * Aul_ / (8 * pi_ * nu_0_ * nu_0_ * (2 * Jl_ + 1));

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
	auto u_dev    = u_ ->view_device();

	auto Doppler_width_dev = Doppler_width_->view_device();
	auto epsilon_dev       = epsilon_ ->view_device();
	auto W_T_dev           = W_T_->view_device();
	auto T_KQ_dev          = T_KQ_->view_device();

	// fill
    Kokkos::parallel_for("INIT-ATM", space_grid_->md_range(), KOKKOS_LAMBDA(int i, int j, int k) 
    {   
    	auto *xi  = xi_dev.block(i, j, k);      
    	auto *T   = T_dev.block(i, j, k);      
    	auto *Nl  = Nl_dev.block(i, j, k);     
    	auto *Cul = Cul_dev.block(i, j, k);      
    	auto *Qel = Qel_dev.block(i, j, k);      
        auto *D2  = D2_dev.block(i, j, k);
        auto *D1  = D1_dev.block(i, j, k);
        auto *k_L = k_L_dev.block(i, j, k);
        auto *u   = u_dev.block(i, j, k);

        auto *Doppler_width = Doppler_width_dev.block(i, j, k);
        auto *epsilon       = epsilon_dev.block(i, j, k);
        auto *W_T           = W_T_dev.block(i, j, k);        

        // precompute quantities depening only on position
        D2[0]  = 0.5 * Qel[0];

		D1[0]  = tmp_const3 * D2[0];

		k_L[0] = tmp_const * Nl[0];

		Doppler_width[0] = dE * std::sqrt(2 * k_B_ * T[0] / mass_real + xi[0] * xi[0]);

		epsilon[0] = Cul[0]/(Cul[0] + Aul_);	

		W_T[0] = tmp_const2 * std::exp(- h_ * nu_0_ / (k_B_ * T[0]));		

		// on position and frequency
		for (size_t n = 0; n < N_nu_; ++n)
		{
			u[n] = (nu_0_ - nu_grid_[n]) / Doppler_width[0];
		}	
    });	

	Kokkos::parallel_for("INIT-TKQ", space_grid_complex_->md_range(), KOKKOS_LAMBDA(int i, int j, int k) 
	{
		auto *T_KQ = T_KQ_dev.block(i, j, k);

		size_t counter = 0;

		// on position and direction
		const int KQ_block_size = 6 * N_dirs_;

		for (int i = 0; i < 4; ++i) // TODO change order?
		{
			for (size_t j = 0; j < N_theta_; ++j)
			{
				for (size_t k = 0; k < N_chi_; ++k)
				{
					auto T_KQ_vec = compute_T_KQ(i, theta_grid_[j], chi_grid_[k]);

					for (int KQ = 0; KQ < 6; ++KQ)
					{
						T_KQ[counter] = T_KQ_vec[KQ];	

						counter++;
					}				
				}
			}
		}
	});	
	   	
	// // precompute etas and rhos //////////////////////////////////////////
	// PetscErrorCode ierr;

	// int istart, iend;
	// int ix[4];

	// bool dichroism_warning = false;

	// ierr = VecGetOwnershipRange(eta_field_, &istart, &iend);CHKERRV(ierr);
	
	// std::vector<Real> etas_and_rhos;

	// size_t i_loc, j_loc = 0, k_loc = 0, n_loc = 0;
	
	// for (int i = istart; i < iend; i = i + 4) 
	// {	
	// 	i_loc = i/block_size_;
		
	// 	etas_and_rhos = compute_eta_and_rhos(i_loc, j_loc, k_loc, n_loc);	

	// 	n_loc++;

	// 	if (n_loc == N_nu_)
	// 	{
	// 		n_loc = 0;
	// 		k_loc++;

	// 		if (k_loc == N_chi_)
	// 		{
	// 			k_loc = 0;
	// 			j_loc++;

	// 			if (j_loc == N_theta_) j_loc = 0;				
	// 		}
	// 	}		

	// 	// assign ix[ii] = i + ii;
	// 	std::iota(ix, ix + 4, i);

	// 	// debug		
	// 	if (etas_and_rhos[0] <= 0)  std::cerr << "\nWARNING: eta_I not positive!" << std::endl; 

	// 	const Real dichroism_module = std::sqrt(etas_and_rhos[1] * etas_and_rhos[1] + etas_and_rhos[2] * etas_and_rhos[2] + etas_and_rhos[3] * etas_and_rhos[3]);
		
	// 	if (etas_and_rhos[0] < dichroism_module) dichroism_warning = true;
	// 	// {			
	// 	// 	std::cout << "eta_I = " << etas_and_rhos[0] << std::endl;
	// 	// 	std::cout << "dichroism_module = " << dichroism_module << std::endl;
	// 	// } 
				
	// 	ierr = VecSetValues(eta_field_, 4, ix, &etas_and_rhos[0], INSERT_VALUES);CHKERRV(ierr); 
	// 	ierr = VecSetValues(rho_field_, 4, ix, &etas_and_rhos[4], INSERT_VALUES);CHKERRV(ierr); 		
	// }

	// ierr = VecAssemblyBegin(eta_field_);CHKERRV(ierr); 
	// ierr = VecAssemblyEnd(eta_field_);CHKERRV(ierr); 
	// ierr = VecAssemblyBegin(rho_field_);CHKERRV(ierr); 
	// ierr = VecAssemblyEnd(rho_field_);CHKERRV(ierr); 

	// if (dichroism_warning) std::cout << "\nWARNING: eta_I < eta! (Eq. (7) Gioele Paganini 2018, Part III)" << std::endl; 	

	// // precomputing dtau /////////////////////////////////////////////////////////////

	// const int stoke_block = block_size_ / 4; // stoke_block = N_chi_ * N_theta * N_nu

	// int eta_start, eta_end;

	// ierr = VecGetOwnershipRange(eta_field_, &eta_start, &eta_end);CHKERRV(ierr);

	// // copy the portion of eta_field_ that is needed //////////
	// const bool not_last_rank = mpi_rank_ < mpi_size_ - 1;

	// // create vector for the first portion of eta own by the next processor
	// Vec next_eta; 		

	// ierr = VecCreate(PETSC_COMM_WORLD, &next_eta);CHKERRV(ierr);	
	// ierr = VecSetSizes(next_eta, stoke_block, mpi_size_ * stoke_block);CHKERRV(ierr);		
	// ierr = VecSetFromOptions(next_eta);CHKERRV(ierr);				
					
	// IS is_send, is_recive;

	// if (not_last_rank)	
	// {
	// 	ierr = ISCreateStride(PETSC_COMM_SELF, stoke_block, eta_end, 4, &is_send);CHKERRV(ierr);
	// 	ierr = ISCreateStride(PETSC_COMM_SELF, stoke_block, stoke_block * mpi_rank_, 1, &is_recive);CHKERRV(ierr);
	// }
	// else
	// {
	// 	ierr = ISCreateStride(PETSC_COMM_SELF, 0, 0, 1, &is_send);CHKERRV(ierr);
	// 	ierr = ISCreateStride(PETSC_COMM_SELF, 0, 0, 1, &is_recive);CHKERRV(ierr);
	// }
		
	// VecScatter scatter_ctx;	

	// ierr = VecScatterCreate(eta_field_,is_send,next_eta,is_recive, &scatter_ctx);CHKERRV(ierr);
	// ierr = VecScatterBegin(scatter_ctx, eta_field_, next_eta, INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
 //    ierr = VecScatterEnd(scatter_ctx, eta_field_, next_eta, INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
    
 //    // clean
 //    ierr = VecScatterDestroy(&scatter_ctx);CHKERRV(ierr);
 //    ierr = ISDestroy(&is_send);CHKERRV(ierr);
	// ierr = ISDestroy(&is_recive);CHKERRV(ierr);

 //    /////////////// scattering done ///////////////
	
 //    ierr = VecGetOwnershipRange(dtau_, &istart, &iend);CHKERRV(ierr);    

 //    int start_next, end_next, counter = 0;
 //    ierr = VecGetOwnershipRange(next_eta, &start_next, &end_next);CHKERRV(ierr);    
          
	// Real eta1, eta2, dtau;
	// int index;		

	// bool sign_err = false; 
	
	// for (int i = istart; i < iend; ++i) 
	// {
	// 	i_loc = i/stoke_block;	
		
	// 	index = i * 4;
	// 	ierr = VecGetValues(eta_field_, 1, &index, &eta1);CHKERRV(ierr);

	// 	index += block_size_;

	// 	if (index < eta_end)
	// 	{
	// 		ierr = VecGetValues(eta_field_, 1, &index, &eta2);	CHKERRV(ierr);

	// 		dtau = 0.5 * (eta1 + eta2) * 1e5 * std::abs(depth_grid_[i_loc + 1] - depth_grid_[i_loc]); // conversion to cm			

	// 		ierr = VecSetValue(dtau_, i, dtau, INSERT_VALUES);CHKERRV(ierr);

	// 		if (dtau <= 0) sign_err = true; 

	// 	}
	// 	else if (not_last_rank)// use data copied from next processor
	// 	{							
	// 		index = start_next + counter;
	// 		ierr = VecGetValues(next_eta, 1, &index, &eta2);CHKERRV(ierr);	
	// 		counter++;

	// 		dtau = 0.5 * (eta1 + eta2) * 1e5 * std::abs(depth_grid_[i_loc + 1] - depth_grid_[i_loc]); // conversion to cm					

	// 		ierr = VecSetValue(dtau_, i, dtau, INSERT_VALUES);CHKERRV(ierr);

	// 		if (dtau <= 0) sign_err = true; 
	// 	}				
	// }

	// if (sign_err) std::cerr << "\nWARNING: dtau negative!" << std::endl; 

	// ierr = VecAssemblyBegin(dtau_);CHKERRV(ierr); 
	// ierr = VecAssemblyEnd(dtau_);CHKERRV(ierr); 

	// // save_vec(dtau_, "../output/dtau.m","dt");

	// // clean	
	// ierr = VecDestroy(&next_eta);CHKERRV(ierr); 		

	if (mpi_rank_ == 0) std::cout << "done" << std::endl;	
}






