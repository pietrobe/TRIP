#include "RT_problem.hpp"


void RT_problem::init_fields(){

		// init fields 
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

void RT_problem::init_atmosphere(){

	// init atmospheric quantities 
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
	u_        = std::make_shared<Field_t>("u", space_grid_, N_nu_, sgrid::BOX_STENCIL);  
	sigma_    = std::make_shared<Field_t>("sigma", space_grid_, N_nu_, sgrid::BOX_STENCIL);
	k_c_      = std::make_shared<Field_t>("k_c", space_grid_, N_nu_, sgrid::BOX_STENCIL);
	eps_c_th_ = std::make_shared<Field_t>("eps_c_th", space_grid_, N_nu_, sgrid::BOX_STENCIL);

	// quantities depending on position and direction
	T_KQ_ = std::make_shared<Field_t>("Doppler_width", space_grid_, N_dirs_, sgrid::BOX_STENCIL);  

	// allocate
	Nl_  -> allocate_on_device(); 
	Nu_  -> allocate_on_device(); 
	T_   -> allocate_on_device(); 
	xi_  -> allocate_on_device(); 
	nu_L_-> allocate_on_device(); 
	Cul_ -> allocate_on_device(); 
	Qel_ -> allocate_on_device(); 
	a_   -> allocate_on_device(); 
	W_T_ -> allocate_on_device(); 
	
	theta_B_ -> allocate_on_device(); 
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




