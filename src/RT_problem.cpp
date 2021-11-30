#include "RT_problem.hpp"


void RT_problem::init_fields(){

		// init fields 
		I_field_ = std::make_shared<Field_t>("I", space_grid_, block_size_, sgrid::BOX_STENCIL);
		S_field_ = std::make_shared<Field_t>("S", space_grid_, block_size_, sgrid::BOX_STENCIL);

		eta_field_ = std::make_shared<Field_t>("eta", space_grid_, block_size_, sgrid::BOX_STENCIL);
		rho_field_ = std::make_shared<Field_t>("rho", space_grid_, block_size_, sgrid::BOX_STENCIL);

		// TODO?
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




