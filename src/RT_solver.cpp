#include "RT_solver.hpp"


void MF_context::apply_bc(Field_ptr_t I_field, const Real I0){

	if (mpi_rank_ == 0) std::cout << "\nApplying BC..." << std::endl;

	// init some quantities 
	const auto N_z        = RT_problem_->N_z_;
	const auto block_size = RT_problem_->block_size_;
	const auto space_grid = RT_problem_->space_grid_;	

	// apply BC
	auto W_T_dev     = RT_problem_->W_T_->view_device();
	auto I_field_dev =           I_field->view_device();
	auto g_dev       =        space_grid->view_device();

	sgrid::parallel_for("APPLY BC", space_grid->md_range(), SGRID_LAMBDA(int i, int j, int k) {
	
		const size_t k_global = g_dev.global_coord(2, k);					

		Real W_T_deep; 

		// just in max depth
		if (k_global == N_z - 1) 
		{		
			W_T_deep = I0 * W_T_dev.ref(i,j,k);

			for (int b = 0; b < (int)block_size; b = b + 4) 
			{
				I_field_dev.block(i,j,k)[b] = W_T_deep; 
			}
		}
	});											
}


void MF_context::find_intersection(double theta, double chi, const double Z_down, const double Z_top, const double L, t_intersect *T) {
    
    // check theta
    if (theta == 0 or chi == 0) std::cout << "WARNING: ray direction not supported!" << std::endl;       
    {        
        theta += 1e-16;
        chi   += 1e-16;                 
    }
    
    // check widths
    if (Z_down <= 0 ) std::cout << "WARNING: Z_down not positive" << std::endl;       
    if (Z_top  <= 0 ) std::cout << "WARNING: Z_top not positive"  << std::endl;       
    if (L      <= 0 ) std::cout << "WARNING: L not positive"      << std::endl;       
    
    // unit vector in the direction of the ray (minus for different convection in formal solver)
    const double x = - sin(theta) * cos(chi);
    const double y = - sin(theta) * sin(chi); 
    const double z = - cos(theta); 

    const int dix = (x > 0. ? 1 : -1);
    const int diy = (y > 0. ? 1 : -1);
    const int diz = (z > 0. ? 1 : -1);

    const double x_inters = L * dix;
    const double y_inters = L * diy;
    const double z_inters = diz > 0 ? Z_top : -Z_down;

    const double tx = x_inters / x;
    const double ty = y_inters / y;
    const double tz = z_inters / z;

    double u1, u2, v1, v2, u, v;

    if (tz <= tx && tz <= ty) {
        T->distance = tz;
        for (int i=0; i<4; i++) T->iz[i] = diz;
        T->ix[0] = MIN(0,dix); T->iy[0] = MIN(0,diy);
        T->ix[1] = MAX(0,dix); T->iy[1] = MIN(0,diy);
        T->ix[2] = MAX(0,dix); T->iy[2] = MAX(0,diy);
        T->ix[3] = MIN(0,dix); T->iy[3] = MAX(0,diy);
        u1 = L * T->ix[0]; u2 = L * T->ix[1];
        v1 = L * T->iy[0]; v2 = L * T->iy[2];
        u = T->distance * x; v = T->distance * y;
    }
    else if (ty <= tx && ty <= tz) {
        T->distance = ty;
        for (int i=0; i<4; i++) T->iy[i] = diy;
        T->ix[0] = MIN(0,dix); T->iz[0] = MIN(0,diz);
        T->ix[1] = MAX(0,dix); T->iz[1] = MIN(0,diz);
        T->ix[2] = MAX(0,dix); T->iz[2] = MAX(0,diz);
        T->ix[3] = MIN(0,dix); T->iz[3] = MAX(0,diz);
        u1 = L * T->ix[0]; u2 = L * T->ix[1];
        v1 = MIN(0., z_inters); v2 = MAX(0., z_inters);
        u = T->distance * x; v = T->distance * z;
    }
    else {
        T->distance = tx;
        for (int i=0; i<4; i++) T->ix[i] = dix;
        T->iy[0] = MIN(0,diy); T->iz[0] = MIN(0,diz);
        T->iy[1] = MAX(0,diy); T->iz[1] = MIN(0,diz);
        T->iy[2] = MAX(0,diy); T->iz[2] = MAX(0,diz);
        T->iy[3] = MIN(0,diy); T->iz[3] = MAX(0,diz);
        u1 = L * T->iy[0]; u2 = L * T->iy[1];
        v1 = MIN(0., z_inters); v2 = MAX(0., z_inters);
        u = T->distance * y; v = T->distance * z;
    }

    double norm = 1.0 / ((u2-u1) * (v2-v1));

    T->w[0] = norm * (u2-u) * (v2-v);
    T->w[1] = norm * (u-u1) * (v2-v);
    T->w[2] = norm * (u-u1) * (v-v1);
    T->w[3] = norm * (u2-u) * (v-v1);

    for (int i = 0; i < 4; ++i)
    {
    	if (T->w[i] < 0 or T->w[i] > 1)
    	{
    		std::cout << "WARNING: w has a problem!" << std::endl;         	
    		std::cout << "theta = " << theta << std::endl;         	
    		std::cout << "chi = "  << chi << std::endl;         	
    		std::cout << "Z_down = "  << Z_down << std::endl;         	
    		std::cout << "Z_top = "  << Z_top << std::endl;         	
    		std::cout << "L = "  << L << std::endl;         	
    		
    		std::cout << "i = " << i << std::endl;         	
    		std::cout << "w = " << T->w[i] << std::endl;         	
    	}  
    }
}


void MF_context::formal_solve_local(Field_ptr_t I_field, const Field_ptr_t S_field, const Real I0){

	 // S_field->exchange_halos(); to include in update S 			

	if (mpi_rank_ == 0) std::cout << "\nStart formal solution..." << std::endl;

	// init some quantities 
	const auto N_theta = RT_problem_->N_theta_;
	const auto N_chi   = RT_problem_->N_chi_;
	const auto N_nu    = RT_problem_->N_nu_;
	const auto N_x     = RT_problem_->N_x_;
	const auto N_y     = RT_problem_->N_y_;
	const auto N_z     = RT_problem_->N_z_;
	const auto N_s     = RT_problem_->N_s_;

	const auto vertical_decompostion = RT_problem_->only_vertical_decompostion_;

	const auto block_size = RT_problem_->block_size_;
	const auto tot_size   = RT_problem_->tot_size_;
	
	const auto mu_grid    = RT_problem_->mu_grid_;
	const auto theta_grid = RT_problem_->theta_grid_;	
	const auto chi_grid   = RT_problem_->chi_grid_;	
	const auto depth_grid = RT_problem_->depth_grid_;	

	const auto L = RT_problem_->L_;	
	const auto space_grid = RT_problem_->space_grid_;	

	const auto eta_dev = RT_problem_->eta_field_->view_device();
	const auto rho_dev = RT_problem_->rho_field_->view_device();

	auto I_dev = 	I_field->view_device();		
	auto S_dev = 	S_field->view_device();	
	auto g_dev = space_grid->view_device();

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	if (i_start > 1 or j_start > 1 or k_start > 1) cout << "WARNING: tested only for margin = 1!"<< endl;

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
	const int k_end = k_start + g_dev.dim[2];	

	int i_aux, j_aux, k_aux, k_global, i_intersect, j_intersect, k_intersect, b_start, b_index;

	// misc coeffs
	double theta, chi, mu, weight, dtau, eta_I_1, Z_down, Z_top;	

	bool boundary;

	// quantities depending on spatial point i
	std::vector<double> I1(4), I2(4), S1(4), S2(4), etas(4), rhos(4), K1(16), K2(16);

	// intersection object
	t_intersect intersection_data;

	// minus for optical depth conversion, trap rule and conversion to cm
	const double coeff = - 0.5 * 1e5;

	for (int rank = 0; rank < mpi_size_; ++rank)
	{			
		// loop over spatial points
		for (int k = k_start; k < k_end; ++k)					
		{													
			// exchange boundary info ---------------------> TODO use buffer to speed up?
			if (not vertical_decompostion) I_field->exchange_halos(); 			    	
			
			for (int j = j_start; j < j_end; ++j)
			{
				for (int i = i_start; i < i_end; ++i)
				{
					// loop over directions (TODO could be parallel)
					for (int j_theta = 0; j_theta < (int)N_theta; ++j_theta)
					{
						theta = theta_grid[j_theta];
						mu    = mu_grid[j_theta];						

						k_aux = (mu > 0.0) ? k_end - k : k; 

						// depth index
						k_global = g_dev.global_coord(2, k_aux);	

						boundary = (k_global == 0 and mu < 0) or (k_global == (int)N_z - 1 and mu > 0);
						
						if (not boundary)
						{						
							// set vertical box size
							Z_down = (mu > 0) ? depth_grid[k_global] -  depth_grid[k_global + 1] : 1.0;
							Z_top  = (mu < 0) ? depth_grid[k_global - 1] - depth_grid[k_global]  : 1.0;						

							for (int k_chi = 0; k_chi < (int)N_chi; ++k_chi)
							{			
								chi = chi_grid[k_chi];
							
								i_aux = (cos(chi) < 0.0) ? i_end - i : i;	
								j_aux = (sin(chi) < 0.0) ? j_end - j : j;

								// loop on freqs
								for (int n = 0; n < (int)N_nu; ++n)
								{
									// block index
									b_start = RT_problem_->local_to_block(j_theta, k_chi, n); 

									// solve ODE
									
									// set S2
									for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
									{				
										b_index = b_start + i_stokes;

										// get eta and rho
										etas[i_stokes] = eta_dev.block(i_aux,j_aux,k_aux)[b_index];						
										rhos[i_stokes] = rho_dev.block(i_aux,j_aux,k_aux)[b_index];		
										
										// set S2 with values in the current grid nodes 							
										S2[i_stokes] = S_dev.block(i_aux,j_aux,k_aux)[b_index];									
									}

									// for optical depth conversion
									eta_I_1 = etas[0];		
								
									K2 = assemble_propagation_matrix_scaled(etas, rhos);
									
									// compute interesection point
									find_intersection(theta, chi, Z_down, Z_top, L, &intersection_data);	

									// set etas, rhos and S1 and I1 to zero
									for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
									{
										etas[i_stokes] = 0;
										rhos[i_stokes] = 0;
										S1[i_stokes]   = 0;
										I1[i_stokes]   = 0;
									}

									// loop over the four vertex of the intersection face
									for (int face_vertices = 0; face_vertices < 4; ++face_vertices)
									{
										i_intersect = i_aux + intersection_data.ix[face_vertices];
										j_intersect = j_aux + intersection_data.iy[face_vertices];
										k_intersect = k_aux - intersection_data.iz[face_vertices]; // minus because k increases going downwards		
																				
										// correct for periodic boundary
										if (vertical_decompostion)
										{
											if (i_intersect == 0)            i_intersect = N_x; 
											if (j_intersect == 0)            j_intersect = N_y;
											if (i_intersect == (int)N_x + 1) i_intersect = 1;
											if (j_intersect == (int)N_y + 1) j_intersect = 1;	
										}										

										weight = intersection_data.w[face_vertices];
										
										for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
										{
											b_index = b_start + i_stokes;										

											// get eta and rho
											etas[i_stokes] += weight * eta_dev.block(i_intersect,j_intersect,k_intersect)[b_index];						
											rhos[i_stokes] += weight * rho_dev.block(i_intersect,j_intersect,k_intersect)[b_index];

											// set S1 and I1
											S1[i_stokes] += weight * S_dev.block(i_intersect,j_intersect,k_intersect)[b_index];												
											I1[i_stokes] += weight * I_dev.block(i_intersect,j_intersect,k_intersect)[b_index];																			
										}
									}								
																	
									K1 = assemble_propagation_matrix_scaled(etas, rhos);
									
									// optical depth step								
									dtau = coeff * (eta_I_1 + etas[0]) * intersection_data.distance;

									if (dtau > 0 ) std::cout << "ERROR in dtau sign" << std::endl;							

									formal_solver_.one_step(dtau, K1, K2, S1, S2, I1, I2);	
									
									// test
									if (j_theta == 1 and k_chi == 0 and n == 0 and rank == mpi_size_ - 1)
									{									
										// std::cout << "I1 = " << I1[0] << std::endl;		
										if (i == 1 and j == 1)
										{
											// std::cout << intersection_data.distance << std::endl;	
											// std::cout << "coeff = " << coeff << std::endl;	
											// std::cout << "eta_I_1 = " << eta_I_1 << std::endl;	
											// std::cout << "etas[0] = " << etas[0] << std::endl;	
											// std::cout << "intersection_data.distance = " << intersection_data.distance << std::endl;					
											
											// std::cout << "I1 = " << I1[0] << std::endl;	
											std::cout << "I2 = " << I2[0] << std::endl;		
											// std::cout << "Z_down = " << Z_down << std::endl;					
											// std::cout << "Z_top = "  << Z_top  << std::endl;					
											// std::cout << "L = "  << L  << std::endl;					
											// std::cout << "theta = "  << theta  << std::endl;	
											// std::cout << "chi = "  << chi  << std::endl;	
											
											// std::cout << "S1 = " << std::endl;
											// for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S1[i_stokes] << std::endl;

											// std::cout << "S2 = " << std::endl;
											// for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S2[i_stokes] << std::endl;
																		
											// std::cout << "K1 = " << std::endl;
											// for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K1[i_stokes] << std::endl;

											// std::cout << "K2 = " << std::endl;
											// for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K2[i_stokes] << std::endl;												
										} 									
									}

									// write result
									for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
									{							
										I_dev.block(i_aux,j_aux,k_aux)[b_start + i_stokes] = I2[i_stokes];										
									}
								}						
							}
						}
					}
				}
			}		
		}
		
		// communication
		if (mpi_size_ > 1) I_field->exchange_halos(); 			
	}
}
 				
// TOOD invert loops order

// void MF_context::make_intensity_positive(Vec &I_field){

// 	PetscErrorCode ierr;

// 	int istart, iend;
// 	ierr = VecGetOwnershipRange(I_field, &istart, &iend);CHKERRV(ierr);

// 	double value;

// 	for (int i = istart; i < iend; i = i + 4)
// 	{
// 		ierr = VecGetValues(I_field, 1, &i, &value);CHKERRV(ierr);	

// 		if (value < 0) ierr = VecSetValue(I_field, i, 0.0, INSERT_VALUES);CHKERRV(ierr);		
// 	}
// }


// void MF_context::set_up_emission_module(){

// 	// Build module
//     auto fsf_in_sh_ptr =
//     rii_include::formal_solver_factory_from_in_struct::make_formal_solver_factory_from_in_struct_shared_ptr();

//     // interface
//     rii_include::in_RT_problem_interface<RT_problem> RT_interface;
//     RT_interface.add_models(RT_problem_, fsf_in_sh_ptr);

//     ecc_sh_ptr_ = rii_include::emission_coefficient_computation::make_emission_coefficient_computation_shared_ptr();

//     bool flag = ecc_sh_ptr_->build_problem(fsf_in_sh_ptr);

//     if (not flag) std::cerr << "Error in set_up_emission_module()!" << std::endl;

//     epsilon_computation_function_ = ecc_sh_ptr_->make_computation_function(
// 		{
// 			rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_R_II,
// 		 	rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_R_III,
// 		 	rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_csc
// 		},  
// 		rii_consts::rii_units::kilometer);

//     epsilon_computation_function_approx_ = ecc_sh_ptr_->make_computation_function(
// 	{
// 		// rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_R_II,
// 	 	rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_R_III_CRD_limit,
// 	 	// rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_csc
// 	},  
// 	rii_consts::rii_units::kilometer);
		 
//     offset_f_ = rii_include::make_default_offset_function(RT_problem_->N_theta_, RT_problem_->N_chi_, RT_problem_->N_nu_);
// }


// // emissivity from scattering (line + continuum)
// void MF_context::update_emission(const Vec &I_field, const bool approx){ 	
	
// 	PetscErrorCode ierr; 
		
// 	// test Simone module
//     const auto block_size = RT_problem_->block_size_;     

//     std::vector<double>  input(block_size);        
//     std::vector<double> output(block_size); 

//     double height;
//     int ix[block_size];

//     int istart, iend; 
//     ierr = VecGetOwnershipRange(I_field, &istart, &iend);CHKERRV(ierr);	
	
// 	const int istart_local = istart / block_size;
// 	const int iend_local   = iend   / block_size;

//     for (int i = istart_local; i < iend_local; ++i)
//     {
//     	// set indeces
//     	std::iota(ix, ix + block_size, i * block_size);

//     	// call Simone module
//     	height = RT_problem_->depth_grid_[i];

//     	// get I field at ith height
//     	ierr = VecGetValues(I_field, block_size, ix, &input[0]);CHKERRV(ierr);	

//     	// set input field
// 		ecc_sh_ptr_->update_incoming_field<double>(height, offset_f_, input.data());

//     	ecc_sh_ptr_->set_threads_number(1);
    	
//     	if (approx)
//     	{
//     		auto epsilon_grid = epsilon_computation_function_approx_(height);
// 	    	rii_include::make_indices_convertion_function<double>(epsilon_grid, offset_f_)(output.data());    
//     	}
//     	else
//     	{
//     		auto epsilon_grid = epsilon_computation_function_(height);
// 	    	rii_include::make_indices_convertion_function<double>(epsilon_grid, offset_f_)(output.data());    
//     	}
    	
//     	// set S_field_ accordingly        	
//     	ierr = VecSetValues(RT_problem_->S_field_, block_size, ix, &output[0], INSERT_VALUES);CHKERRV(ierr);		
//     }
    
// 	// test S = I;
// 	// ierr = VecCopy(I_field, RT_problem_->S_field_);CHKERRV(ierr);		

// 	// test S = 1
// 	// ierr = VecSet(RT_problem_->S_field_, 1.0);CHKERRV(ierr);		

// 	// test S = 0;
// 	// ierr = VecZeroEntries(RT_problem_->S_field_);CHKERRV(ierr);		

// 	// // LTE S = W_T;
// 	// const auto block_size = RT_problem_->block_size_;
 		
// 	// int istart, iend; 
// 	// ierr = VecGetOwnershipRange(RT_problem_->I_field_, &istart, &iend);CHKERRV(ierr);	

// 	// size_t i_space;
// 	// double value;

// 	// for (int i = istart; i < iend; ++i)
// 	// {
// 	// 	i_space = i / block_size;

// 	// 	value = RT_problem_->W_T_[i_space];

// 	// 	ierr = VecSetValues(RT_problem_->S_field_, 1, &i, &value, INSERT_VALUES);CHKERRV(ierr); 							
// 	// }

// 	// ierr = VecAssemblyBegin(RT_problem_->S_field_);CHKERRV(ierr); 
//  	// ierr = VecAssemblyEnd(RT_problem_->S_field_);CHKERRV(ierr); 
    
//  	// scale emission by eta_I 
// 	int ixx[4];
// 	double eta_I, eta_I_inv;	
// 	double S[4];
		
// 	for (int i = istart; i < iend; i = i + 4)
// 	{
// 		ierr = VecGetValues(RT_problem_->eta_field_, 1, &i, &eta_I);CHKERRV(ierr);	

// 		eta_I_inv = 1.0 / eta_I;

// 		std::iota(ixx, ixx + 4, i);

// 		ierr = VecGetValues(RT_problem_->S_field_, 4, ixx, S);CHKERRV(ierr);	

// 		S[0] = eta_I_inv * S[0];
// 		S[1] = eta_I_inv * S[1];
// 		S[2] = eta_I_inv * S[2];
// 		S[3] = eta_I_inv * S[3];

// 		ierr = VecSetValues(RT_problem_->S_field_, 4, ixx, S, INSERT_VALUES);CHKERRV(ierr); 									
// 	}

// 	ierr = VecAssemblyBegin(RT_problem_->S_field_);CHKERRV(ierr); 
//  	ierr = VecAssemblyEnd(RT_problem_->S_field_);CHKERRV(ierr);  
// }

// void RT_solver::save_Lamda(){

// 	if (mpi_size_ > 1) std::cout << "WARNING in save mat!"<< std::endl;

// 	const bool L_flag = false;

// 	PetscErrorCode ierr;

// 	const auto tot_size = RT_problem_->tot_size_;

// 	Mat L;

// 	ierr = MatCreate(PETSC_COMM_WORLD, &L);CHKERRV(ierr);
//     ierr = MatSetSizes(L,PETSC_DECIDE, PETSC_DECIDE, tot_size, tot_size);CHKERRV(ierr);
//     ierr = MatSetType(L, MATDENSE     );CHKERRV(ierr);
//     ierr = MatSetUp(L);CHKERRV(ierr);

// 	Vec d_i;

// 	ierr = VecCreate(PETSC_COMM_WORLD, &d_i);CHKERRV(ierr);
//     ierr = VecSetSizes(d_i,PETSC_DECIDE, tot_size);CHKERRV(ierr);    
//     ierr = VecSetType(d_i, VECSEQ);CHKERRV(ierr);    

//     int ix[tot_size];
//   	std::iota(ix, ix + tot_size, 0);

//   	double * vals[tot_size];  	
//   	if (L_flag)
//   	{
//   		ierr = VecGetArray(RT_problem_->I_field_, vals);CHKERRV(ierr);  		
//   	}
//   	else
//   	{
//   		ierr = VecGetArray(RT_problem_->S_field_, vals);CHKERRV(ierr);
//   	}
  	
// 	// fill L rows
// 	for (int i = 0; i < tot_size; ++i)
// 	{

// 		std::cout << i << std::endl;

// 		ierr = VecZeroEntries(d_i);CHKERRV(ierr);               
// 		ierr = VecSetValue(d_i, i, 1.0, INSERT_VALUES);CHKERRV(ierr);
	    
// 		if (L_flag)
// 	  	{
// 	  		mf_ctx_.formal_solve(RT_problem_->I_field_, d_i, 0.0);
// 	  	}
// 	  	else
// 	  	{
// 	  		mf_ctx_.update_emission(d_i);  
// 	  	}
  	    	
// 	  	ierr = MatSetValues(L,1,&i,tot_size,ix,*vals,INSERT_VALUES);CHKERRV(ierr);
// 	}

// 	ierr = MatAssemblyBegin(L,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
//     ierr = MatAssemblyEnd(L,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);

// 	if (L_flag)
//   	{
//   		save_mat(L, "../output/L.m" ,"Lam"); 
//   	}
//   	else
//   	{
//   		save_mat(L, "../output/S.m" ,"Sigma"); 
//   	}

// 	ierr = MatDestroy(&L);CHKERRV(ierr);
// 	ierr = VecDestroy(&d_i);CHKERRV(ierr);
// };

// void RT_solver::assemble_rhs(){

//   	if (mpi_rank_ == 0) std::cout << "\nAssembling right hand side...";

// 	PetscErrorCode ierr;

// 	const auto N_nu       = RT_problem_->N_nu_;
// 	const auto       size = RT_problem_->tot_size_;
// 	const auto block_size = RT_problem_->block_size_;
// 	const auto eta_field  = RT_problem_->eta_field_;

// 	Vec eps_th; // TODO use directly rhs_
// 	double eta_i, eta_I, value;
// 	size_t i_space;
// 	int index_I, index_s_nu;

// 	ierr = VecCreate(PETSC_COMM_WORLD, &eps_th);CHKERRV(ierr);    
// 	ierr = VecSetSizes(eps_th,PETSC_DECIDE,size);CHKERRV(ierr);   
// 	ierr = VecSetBlockSize(eps_th,block_size);CHKERRV(ierr);
// 	ierr = VecSetFromOptions(eps_th);CHKERRV(ierr);

// 	int istart, iend;
// 	ierr = VecGetOwnershipRange(eta_field, &istart, &iend);CHKERRV(ierr);

// 	// fill eps_th =  eps_c_th +  eps_l_th
// 	for (int i = istart; i < iend; ++i)
// 	{	
// 		value = 0.0;

// 		std::vector<size_t> local_idx = RT_problem_->local_to_global(i);
// 		i_space = local_idx[0];

// 		if (LTE_) // eps = W_T
// 		{
// 			if (local_idx[4] == 0) value = RT_problem_->W_T_[i_space]; // TODO divide by eta_I?
// 		}
// 		else // eps = (eps_c_th + eps_l_th_) / eta_I
// 		{
// 			// get eta_i		
// 			ierr = VecGetValues(eta_field, 1, &i, &eta_i);
			
// 			// first Stokes parameter
// 			if (local_idx[4] == 0)				
// 			{	
// 				index_s_nu = N_nu * i_space + local_idx[3];
				
// 				if (RT_problem_->enable_continuum_) 
// 				{
// 					// eps_c_th
// 					value = (RT_problem_->eps_c_th_[index_s_nu]);	

// 					// eps_l_th		
// 					value += (RT_problem_->epsilon_[i_space]) * (RT_problem_->W_T_[i_space]) * (eta_i - RT_problem_->k_c_[index_s_nu]);				
// 				}
// 				else
// 				{
// 					// eps_l_th
// 					value += (RT_problem_->epsilon_[i_space]) * (RT_problem_->W_T_[i_space]) * eta_i;				
// 				}

// 				// (eps_c_th + eps_l_th_) / eta_I
// 				value /= eta_i;
// 			} 
// 			else
// 			{
// 				// get eta_I (!= eta_i)
// 				index_I = i - local_idx[4];
// 				ierr = VecGetValues(eta_field, 1, &index_I, &eta_I);

// 				// eps_l_th / eta_i_l
// 				value = eta_i * (RT_problem_->epsilon_[i_space]) * (RT_problem_->W_T_[i_space]) / eta_I;				
// 			}			
// 		}
		
// 		ierr = VecSetValue(eps_th,i,value,INSERT_VALUES);CHKERRV(ierr);   
// 	}

// 	ierr = VecAssemblyBegin(eps_th);CHKERRV(ierr); 
// 	ierr = VecAssemblyEnd(eps_th);CHKERRV(ierr);

// 	// save_vec(eps_th, "../output/epsth.m" ,"eps_th"); 
	
// 	// rhs assembly
// 	ierr = VecCreate(PETSC_COMM_WORLD, &rhs_);CHKERRV(ierr);    
// 	ierr = VecSetSizes(rhs_,PETSC_DECIDE,size);CHKERRV(ierr);   
// 	ierr = VecSetBlockSize(rhs_,block_size);CHKERRV(ierr);
// 	ierr = VecSetFromOptions(rhs_);CHKERRV(ierr);

// 	// fill rhs_ from formal solve with bc
// 	mf_ctx_.formal_solve(rhs_, eps_th, 1.0); 	
	
// 	// clean
// 	ierr = VecDestroy(&eps_th);CHKERRV(ierr);

// 	if (mpi_rank_ == 0) std::cout << "done" << std::endl;	
// }


// // set initial guess for first stokes parameter
// void RT_solver::set_I_from_input(const std::string input_path, Vec &I)
// {
// 	if (mpi_rank_ == 0) std::cout << "Reading input initial guess from " << input_path << "/StokesI" << std::endl;

// 	PetscErrorCode ierr;
	
// 	const auto N_theta = RT_problem_->N_theta_;
// 	const auto N_chi   = RT_problem_->N_chi_;
// 	const auto N_nu    = RT_problem_->N_nu_;

// 	const auto block_size = RT_problem_->block_size_;

// 	int istart, iend;
// 	ierr = VecGetOwnershipRange(I, &istart, &iend);CHKERRV(ierr);	

// 	const int istart_local = istart / block_size;
// 	const int iend_local   = iend  / block_size;

// 	std::string filename;
	
// 	for (int i = istart_local; i < iend_local; ++i)
// 	{
// 		std::stringstream ss;

// 		if (i < 10)
// 		{
// 			ss << input_path << "/StokesI/StokesI_00" << i << ".dat";			
// 		}	
// 		else
// 		{
// 			ss << input_path << "/StokesI/StokesI_0" << i << ".dat";
// 		}

// 		filename = ss.str();
		
// 		std::ifstream myFile(filename);

// 		if (not myFile.good()) std::cerr << "\nERROR: File " << filename << " does not exist!\n" << std::endl;

// 		// read each line (size N_nu)
// 		std::string line;	

// 		// direction indeces (for each line)
// 		size_t j, k, direction = 0;

// 		// global position and value;
// 		size_t global_index;
// 		double I_value;			

// 		while(getline(myFile, line))
// 		{
// 			if (direction >= N_theta * N_chi)  std::cerr << "\nERROR input file: direction >= N_theta * N_chi!" << std::endl;

// 			std::istringstream lineStream(line);

// 			for (size_t n = 0; n < N_nu; ++n)
// 			{				
// 				lineStream >> I_value;	

// 				j = direction / N_chi;
// 				k = direction % N_chi;

// 				global_index = RT_problem_->local_to_global(i, j, k, n);

// 				ierr = VecSetValue(I, global_index, I_value, INSERT_VALUES);CHKERRV(ierr);						
// 			}						

// 			direction++;			
// 		} 				
// 	}	

// 	ierr = VecAssemblyBegin(I);CHKERRV(ierr); 
// 	ierr = VecAssemblyEnd(I);CHKERRV(ierr); 
// }

// // matrix-free matrix vector multiplication y = (Id - LJ)x
// PetscErrorCode UserMult(Mat mat,Vec x,Vec y){

// 	PetscErrorCode ierr; 

// 	void *ptr;
//    	MatShellGetContext(mat,&ptr);
//   	MF_context *mf_ctx_ = (MF_context *)ptr;

//   	// compute new emission in S_field_
//   	mf_ctx_->update_emission(x);   

//   	// fill rhs_ from formal solve with zero bc
//   	if (mf_ctx_->threaded_)
//   	{
//   		mf_ctx_->formal_solve_threaded(y, mf_ctx_->RT_problem_->S_field_, 0.0);  		
//   	}
//   	else
//   	{
// 		mf_ctx_->formal_solve(y, mf_ctx_->RT_problem_->S_field_, 0.0);
//   	}

// 	// update I_out = I_in - I_fs (y = x - y)
// 	ierr = VecAYPX(y, -1.0, x);CHKERRQ(ierr);

//   	return ierr;
// }


// // matrix-free matrix vector multiplication y = (Id - LJ)x
// PetscErrorCode UserMult_approx(Mat mat,Vec x,Vec y){

// 	PetscErrorCode ierr; 

// 	void *ptr;
//    	MatShellGetContext(mat,&ptr);
//   	MF_context *mf_ctx_ = (MF_context *)ptr;

//   	// compute new emission in S_field_
//   	mf_ctx_->update_emission(x, true); 
  	
//   	// fill rhs_ from formal solve with zero bc
// 	mf_ctx_->formal_solve(y, mf_ctx_->RT_problem_->S_field_, 0.0);
	
// 	// update I_out = I_in - I_fs (y = x - y)
// 	ierr = VecAYPX(y, -1.0, x);CHKERRQ(ierr);

//   	return ierr;
// }


// PetscErrorCode MF_pc_Destroy(PC pc){

// 	PetscErrorCode ierr;

// 	MF_context *mf_ctx;

// 	ierr = PCShellGetContext(pc,(void**)&mf_ctx); CHKERRQ(ierr);   
//     ierr = PetscFree(mf_ctx);

//     // TODO destroy?

// 	return ierr;

// }

// PetscErrorCode MF_pc_Apply(PC pc,Vec x,Vec y){

// 	PetscErrorCode ierr;

// 	MF_context *mf_ctx;

// 	ierr = PCShellGetContext(pc,(void**)&mf_ctx);CHKERRQ(ierr);   

// 	// apply	
// 	ierr = KSPSolve(mf_ctx->pc_solver_, x, y);CHKERRQ(ierr);

// 	return ierr;
// }
