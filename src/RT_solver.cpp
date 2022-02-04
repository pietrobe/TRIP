#include "RT_solver.hpp"
#include "sgrid_SliceHalo.hpp"


void MF_context::field_to_vec(const Field_ptr_t field, Vec &v)
{
	if (mpi_rank_ == 0) std::cout << "\nCopying field to Vec..." << std::endl;

	PetscErrorCode ierr; 
	
	const auto N_x = RT_problem_->N_x_;
	const auto N_y = RT_problem_->N_y_;
	const auto N_z = RT_problem_->N_z_;

	const auto space_grid = RT_problem_->space_grid_;	
	const auto block_size = RT_problem_->block_size_;

	auto g_dev = space_grid->view_device();
	auto f_dev =      field->view_device();	

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
	const int k_end = k_start + g_dev.dim[2];	

	int istart, iend, row;

	double value;
	
	ierr = VecGetOwnershipRange(v, &istart, &iend);CHKERRV(ierr);	

	int counter = 0;

	for (int k = k_start; k < k_end; ++k)					
	{															
		for (int j = j_start; j < j_end; ++j)
		{
			for (int i = i_start; i < i_end; ++i)				
			{
				for (int b = 0; b < (int)block_size; b = b + 4) 
				{
					// set row index and correposnding entry
					row = istart + counter;

					value = f_dev.block(i, j, k)[b];

					ierr = VecSetValue(v, row, f_dev.block(i,j,k)[b], INSERT_VALUES);CHKERRV(ierr); 

					counter++;
				}							
			}
		}
	}

	ierr = VecAssemblyBegin(v);CHKERRV(ierr); 
	ierr = VecAssemblyEnd(v);CHKERRV(ierr); 
}


void MF_context::vec_to_field(const Field_ptr_t field, Vec &v)
{
	if (mpi_rank_ == 0) std::cout << "\nCopying Vec to field..." << std::endl;

	PetscErrorCode ierr; 
	
	const auto N_x = RT_problem_->N_x_;
	const auto N_y = RT_problem_->N_y_;
	const auto N_z = RT_problem_->N_z_;

	const auto space_grid = RT_problem_->space_grid_;	
	const auto block_size = RT_problem_->block_size_;

	auto g_dev = space_grid->view_device();
	auto f_dev =      field->view_device();	

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
	const int k_end = k_start + g_dev.dim[2];	

	int istart, iend, row;

	double value;
	
	ierr = VecGetOwnershipRange(v, &istart, &iend);CHKERRV(ierr);	

	int counter = 0;

	for (int k = k_start; k < k_end; ++k)					
	{															
		for (int j = j_start; j < j_end; ++j)
		{
			for (int i = i_start; i < i_end; ++i)				
			{
				for (int b = 0; b < (int)block_size; b = b + 4) 
				{
					// set row index and correposnding entry
					row = istart + counter;					

					ierr = VecGetValues(v, 1, &row, &value);CHKERRV(ierr); 

					f_dev.block(i,j,k)[b] = value;								

					counter++;
				}							
			}
		}
	}
}


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


void MF_context::get_2D_weigths(const double x, const double y, double *w)
{
	// get weigths in the unit square
	if (x >= 1 or x <= 0) cout << "Problem in x input "<< endl;
	if (y >= 1 or y <= 0) cout << "Problem in y input "<< endl;

	const double xy = x * y;

	w[0] = 1.0 - x - y + xy; // (1.0 - x) * (1.0 - y) 
	w[1] = x - xy; // (1.0 - y) * x
	w[2] = xy;
	w[3] = y - xy;  //(1.0 - x) * y; 
}

void MF_context::formal_solve_local(Field_ptr_t I_field, const Field_ptr_t S_field, const Real I0){

	// S_field->exchange_halos(); to include in update S 			

	if (mpi_rank_ == 0) std::cout << "\nStart local formal solution..." << std::endl;

	// init some quantities 
	const auto N_theta = RT_problem_->N_theta_;
	const auto N_chi   = RT_problem_->N_chi_;
	const auto N_nu    = RT_problem_->N_nu_;
	const auto N_x     = RT_problem_->N_x_;
	const auto N_y     = RT_problem_->N_y_;
	const auto N_z     = RT_problem_->N_z_;
	const auto N_s     = RT_problem_->N_s_;

	const auto vertical_decomposition = RT_problem_->only_vertical_decomposition_;

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

	const int block_size_half = block_size/2;

	// Initialize slice halo handler
    sgrid::SliceHalo<Field_t> halos_xy(*I_field, 2);   

    // impose boundary conditions
    apply_bc(I_field, I0);     

	for (int rank = 0; rank < mpi_size_; ++rank)
	{			
		// loop over spatial points
		for (int k = k_start; k < k_end; ++k)					
		{													
			// exchange (i,j)-plane boundary info 
			if (not vertical_decomposition) 
			{
				halos_xy.exchange(k);
				halos_xy.exchange(k_end - k); 	
			}			
			
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

								// compute interesection point and init weights
								find_intersection(theta, chi, Z_down, Z_top, L, &intersection_data);	
					
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
																		
									// set etas, rhos and S1 and I1 to zero
									for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
									{									
										if (use_log_interpolation_)
										{
											etas[i_stokes] = 1;
											rhos[i_stokes] = 1;
											S1[i_stokes]   = 1;
											I1[i_stokes]   = 1;											
										}
										else // linear
										{
											etas[i_stokes] = 0;
											rhos[i_stokes] = 0;
											S1[i_stokes]   = 0;
											I1[i_stokes]   = 0;
										}																	
									}

									// loop over the four vertex of the intersection face
									for (int face_vertices = 0; face_vertices < 4; ++face_vertices)
									{
										i_intersect = i_aux + intersection_data.ix[face_vertices];
										j_intersect = j_aux + intersection_data.iy[face_vertices];
										k_intersect = k_aux - intersection_data.iz[face_vertices]; // minus because k increases going downwards		
																				
										// correct for periodic boundary
										if (vertical_decomposition)
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

											if (use_log_interpolation_)
											{
												// get eta and rho
												etas[i_stokes] *= pow(eta_dev.block(i_intersect,j_intersect,k_intersect)[b_index], weight); 
												rhos[i_stokes] *= pow(rho_dev.block(i_intersect,j_intersect,k_intersect)[b_index], weight);

												// set S1 and I1
												S1[i_stokes] *= pow(S_dev.block(i_intersect,j_intersect,k_intersect)[b_index], weight);												
												I1[i_stokes] *= pow(I_dev.block(i_intersect,j_intersect,k_intersect)[b_index], weight);																														
											}
											else
											{
												// get eta and rho
												etas[i_stokes] += weight * eta_dev.block(i_intersect,j_intersect,k_intersect)[b_index]; 
												rhos[i_stokes] += weight * rho_dev.block(i_intersect,j_intersect,k_intersect)[b_index];

												// set S1 and I1
												S1[i_stokes] += weight * S_dev.block(i_intersect,j_intersect,k_intersect)[b_index];												
												I1[i_stokes] += weight * I_dev.block(i_intersect,j_intersect,k_intersect)[b_index];																														
											}									
										}
									}								
																	
									K1 = assemble_propagation_matrix_scaled(etas, rhos);
									
									// optical depth step								
									dtau = coeff * (eta_I_1 + etas[0]) * intersection_data.distance;

									if (dtau > 0 ) std::cout << "ERROR in dtau sign" << std::endl;							

									formal_solver_.one_step(dtau, K1, K2, S1, S2, I1, I2);	
									
									// // test
									// if (j_theta == N_theta/2 and k_chi == 0 and n == 0 and rank == mpi_size_ - 1)
									// {									
									// 	// std::cout << "I1 = " << I1[0] << std::endl;		
									// 	if (i == 1 and j == 1 and k == k_end - 1)
									// 	{
									// 		// std::cout << intersection_data.distance << std::endl;	
									// 		// std::cout << "coeff = " << coeff << std::endl;	
									// 		// std::cout << "eta_I_1 = " << eta_I_1 << std::endl;	
									// 		// std::cout << "etas[0] = " << etas[0] << std::endl;	
									// 		// std::cout << "intersection_data.distance = " << intersection_data.distance << std::endl;					
											
									// 		std::cout << "I1 = " << I1[0] << std::endl;	
									// 		std::cout << "I2 = " << I2[0] << std::endl;		
									// 		// std::cout << "Z_down = " << Z_down << std::endl;					
									// 		// std::cout << "Z_top = "  << Z_top  << std::endl;					
									// 		// std::cout << "L = "  << L  << std::endl;					
									// 		// std::cout << "theta = "  << theta  << std::endl;	
									// 		// std::cout << "chi = "  << chi  << std::endl;	
											
									// 		// std::cout << "S1 = " << std::endl;
									// 		// for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S1[i_stokes] << std::endl;

									// 		// std::cout << "S2 = " << std::endl;
									// 		// for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S2[i_stokes] << std::endl;
																		
									// 		// std::cout << "K1 = " << std::endl;
									// 		// for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K1[i_stokes] << std::endl;

									// 		// std::cout << "K2 = " << std::endl;
									// 		// for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K2[i_stokes] << std::endl;												
									// 	} 									
									// }

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
 				

void MF_context::formal_solve_global(Field_ptr_t I_field, const Field_ptr_t S_field, const Real I0)
{

	// S_field->exchange_halos(); to include in update S 			

	if (mpi_rank_ == 0) std::cout << "\nStart global formal solution..." << std::endl;

	// init some quantities 
	const auto N_theta = RT_problem_->N_theta_;
	const auto N_chi   = RT_problem_->N_chi_;
	const auto N_nu    = RT_problem_->N_nu_;
	const auto N_x     = RT_problem_->N_x_;
	const auto N_y     = RT_problem_->N_y_;
	const auto N_z     = RT_problem_->N_z_;
	const auto N_s     = RT_problem_->N_s_;

	const auto vertical_decomposition = RT_problem_->only_vertical_decomposition_;

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
	double theta, chi, mu, weight, dtau, eta_I_1, Z_down, Z_top, cos_chi, sin_chi;	

	// coeffs for lon ray 
	double phi, R, tmp_coeff2, tmp_coeff3;

	bool boundary, horizontal_face;

	// quantities depending on spatial point i
	std::vector<double> I1(4), I2(4), S1(4), S2(4), etas(4), rhos(4), K1(16), K2(16);

	// intersection object
	t_intersect intersection_data;

	// minus for optical depth conversion, trap rule and conversion to cm
	const double coeff = - 0.5 * 1e5;

	const double pi_half = 0.5 * PI;

	const int block_size_half = 0.5 * block_size;

	// Initialize slice halo handler
    sgrid::SliceHalo<Field_t> halos_xy(*I_field, 2);    

    // impose boundary conditions
    apply_bc(I_field, I0);     
    
	for (int rank = 0; rank < mpi_size_; ++rank)
	{			
		// loop over spatial points
		for (int k = k_start; k < k_end; ++k)					
		{						
			// exchange (i,j)-plane boundary info 
			if (not vertical_decomposition) 
			{
				halos_xy.exchange(k);
				halos_xy.exchange(k_end - k); 	
			}			

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

								cos_chi = cos(chi);
								sin_chi = sin(chi);
							
								i_aux = (cos_chi < 0.0) ? i_end - i : i;	
								j_aux = (sin_chi < 0.0) ? j_end - j : j;		

								find_intersection(theta, chi, Z_down, Z_top, L, &intersection_data);

								horizontal_face = intersection_data.iz[0] == intersection_data.iz[1] and 
									     		  intersection_data.iz[0] == intersection_data.iz[2] and 
									     		  intersection_data.iz[0] == intersection_data.iz[3];

								const bool long_ray = not horizontal_face and (i == i_start or j == j_start);								     		  	
								
								// r = how many blocks are traversed by the current ray								
								double r = 1.0;

								// check if a vertical face is intersected and use long ray
								if (long_ray)																																															
								{																	
									phi = abs(pi_half - theta);

									const double distance = (mu > 0) ? abs(Z_down / sin(phi)) : abs(Z_top / sin(phi));																	
									
									R = distance * cos(phi);

									r = R / L;									
									tmp_coeff2 = cos_chi * r;
									tmp_coeff3 = sin_chi * r;

									// check if long ray goes to an other processor
									const bool i_in_other_proc = (i + floor(tmp_coeff2) < i_start) or (i + ceil( tmp_coeff2) >= i_end);
									const bool j_in_other_proc = (j + floor(tmp_coeff3) < j_start) or (j + ceil( tmp_coeff3) >= j_end);

									if (i_in_other_proc or j_in_other_proc)									
									{
										// nothing: use short ray	
									}
									else // rewrite t_intersect data
									{
										intersection_data.distance = distance;

										intersection_data.ix[0] = floor(tmp_coeff2);
										intersection_data.ix[1] = ceil( tmp_coeff2);
										intersection_data.ix[2] = floor(tmp_coeff2);
										intersection_data.ix[3] = ceil( tmp_coeff2);

										intersection_data.iy[0] = floor(tmp_coeff3);
										intersection_data.iy[1] = floor(tmp_coeff3);
										intersection_data.iy[2] = ceil( tmp_coeff3);
										intersection_data.iy[3] = ceil( tmp_coeff3);										

										for (int face_vertices = 0; face_vertices < 4; ++face_vertices)
										{											
											intersection_data.iz[face_vertices] = (mu > 0) ? -1 : 1;
										
											const double x = abs(tmp_coeff2) - floor(abs(tmp_coeff2));
											const double y = abs(tmp_coeff3) - floor(abs(tmp_coeff3));

											get_2D_weigths(x, y, intersection_data.w);		
										}																			
									}																	
								}									

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
																									
									// set etas, rhos and S1 and I1 to zero									
									for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
									{									
										if (use_log_interpolation_)
										{
											etas[i_stokes] = 1;
											rhos[i_stokes] = 1;
											S1[i_stokes]   = 1;
											I1[i_stokes]   = 1;											
										}
										else // linear
										{
											etas[i_stokes] = 0;
											rhos[i_stokes] = 0;
											S1[i_stokes]   = 0;
											I1[i_stokes]   = 0;
										}																	
									}
									
									// loop over the four vertex of the intersection face
									for (int face_vertices = 0; face_vertices < 4; ++face_vertices)
									{
										i_intersect = i_aux + intersection_data.ix[face_vertices];
										j_intersect = j_aux + intersection_data.iy[face_vertices];
										k_intersect = k_aux - intersection_data.iz[face_vertices]; // minus because k increases going downwards		
																			
										// correct for periodic boundary
										if (vertical_decomposition)
										{
											i_intersect = i_intersect % N_x; 
											j_intersect = j_intersect % N_y;
											if (i_intersect == 0) i_intersect = N_x;
											if (j_intersect == 0) j_intersect = N_y;										
										}	 																		

										weight = intersection_data.w[face_vertices];
									
										for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
										{
											b_index = b_start + i_stokes;										

											if (use_log_interpolation_)
											{
												// get eta and rho
												etas[i_stokes] *= pow(eta_dev.block(i_intersect,j_intersect,k_intersect)[b_index], weight); 
												rhos[i_stokes] *= pow(rho_dev.block(i_intersect,j_intersect,k_intersect)[b_index], weight);

												// set S1 and I1
												S1[i_stokes] *= pow(S_dev.block(i_intersect,j_intersect,k_intersect)[b_index], weight);												
												I1[i_stokes] *= pow(I_dev.block(i_intersect,j_intersect,k_intersect)[b_index], weight);																														
											}
											else
											{
												// get eta and rho
												etas[i_stokes] += weight * eta_dev.block(i_intersect,j_intersect,k_intersect)[b_index]; 
												rhos[i_stokes] += weight * rho_dev.block(i_intersect,j_intersect,k_intersect)[b_index];

												// set S1 and I1
												S1[i_stokes] += weight * S_dev.block(i_intersect,j_intersect,k_intersect)[b_index];												
												I1[i_stokes] += weight * I_dev.block(i_intersect,j_intersect,k_intersect)[b_index];	
											}									
										}											
									}										
																	
									K1 = assemble_propagation_matrix_scaled(etas, rhos);									
									
									// optical depth step								
									dtau = coeff * (eta_I_1 + etas[0]) * intersection_data.distance;									

									if (dtau > 0 ) std::cout << "ERROR in dtau sign" << std::endl;										

									formal_solver_.one_step(dtau, K1, K2, S1, S2, I1, I2);	
									
									// // test
									// if (j_theta == N_theta/2 and k_chi == 0 and n == 0 and rank == mpi_size_ - 1)
									// {									
									// 	// std::cout << "I1 = " << I1[0] << std::endl;		
									// 	if (i == 1 and j == 1 and k == k_end - 1)
									// 	{
									// 		// std::cout << intersection_data.distance << std::endl;	
									// 		// std::cout << "coeff = " << coeff << std::endl;	
									// 		// std::cout << "eta_I_1 = " << eta_I_1 << std::endl;	
									// 		// std::cout << "etas[0] = " << etas[0] << std::endl;	
									// 		// std::cout << "intersection_data.distance = " << intersection_data.distance << std::endl;					
											
									// 		std::cout << "I1 = " << I1[0] << std::endl;	
									// 		std::cout << "I2 = " << I2[0] << std::endl;		
									// 		// std::cout << "Z_down = " << Z_down << std::endl;					
									// 		// std::cout << "Z_top = "  << Z_top  << std::endl;					
									// 		// std::cout << "L = "  << L  << std::endl;					
									// 		// std::cout << "theta = "  << theta  << std::endl;	
									// 		// std::cout << "chi = "  << chi  << std::endl;	
											
									// 		// std::cout << "S1 = " << std::endl;
									// 		// for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S1[i_stokes] << std::endl;

									// 		// std::cout << "S2 = " << std::endl;
									// 		// for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S2[i_stokes] << std::endl;
																		
									// 		// std::cout << "K1 = " << std::endl;
									// 		// for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K1[i_stokes] << std::endl;

									// 		// std::cout << "K2 = " << std::endl;
									// 		// for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K2[i_stokes] << std::endl;												
									// 	} 									
									// }									

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
		if (mpi_size_ > 1) I_field->exchange_halos(); // --------> TODO	better
	}
}

	
void MF_context::set_up_emission_module(){

	// // Build module
 //    auto fsf_in_sh_ptr =
 //    rii_include::formal_solver_factory_from_in_struct::make_formal_solver_factory_from_in_struct_shared_ptr();

 //    // interface
 //    rii_include::in_RT_problem_interface<RT_problem> RT_interface;
 //    RT_interface.add_models(RT_problem_, fsf_in_sh_ptr);

 //    ecc_sh_ptr_ = rii_include::emission_coefficient_computation::make_emission_coefficient_computation_shared_ptr();

 //    bool flag = ecc_sh_ptr_->build_problem(fsf_in_sh_ptr);

 //    if (not flag) std::cerr << "Error in set_up_emission_module()!" << std::endl;

 //    epsilon_computation_function_ = ecc_sh_ptr_->make_computation_function(
	// 	{
	// 		rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_R_II,
	// 	 	rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_R_III,
	// 	 	rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_csc
	// 	},  
	// 	rii_consts::rii_units::kilometer);

 //    epsilon_computation_function_approx_ = ecc_sh_ptr_->make_computation_function(
	// {
	// 	// rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_R_II,
	//  	rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_R_III_CRD_limit,
	//  	// rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_csc
	// },  
	// rii_consts::rii_units::kilometer);
		 
 //    offset_f_ = rii_include::make_default_offset_function(RT_problem_->N_theta_, RT_problem_->N_chi_, RT_problem_->N_nu_);
}


// emissivity from scattering (line + continuum)
void MF_context::update_emission(const Vec &I_field, const bool approx){ 	
	
	PetscErrorCode ierr; 
		
	// // test Simone module
 //    const auto block_size = RT_problem_->block_size_;     

 //    std::vector<double>  input(block_size);        
 //    std::vector<double> output(block_size); 

 //    double height;
 //    int ix[block_size];

 //    int istart, iend; 
 //    ierr = VecGetOwnershipRange(I_field, &istart, &iend);CHKERRV(ierr);	
	
	// const int istart_local = istart / block_size;
	// const int iend_local   = iend   / block_size;

 //    for (int i = istart_local; i < iend_local; ++i)
 //    {
 //    	// set indeces
 //    	std::iota(ix, ix + block_size, i * block_size);

 //    	// call Simone module
 //    	height = RT_problem_->depth_grid_[i];

 //    	// get I field at ith height
 //    	ierr = VecGetValues(I_field, block_size, ix, &input[0]);CHKERRV(ierr);	

 //    	// set input field
	// 	ecc_sh_ptr_->update_incoming_field<double>(height, offset_f_, input.data());

 //    	ecc_sh_ptr_->set_threads_number(1);
    	
 //    	if (approx)
 //    	{
 //    		auto epsilon_grid = epsilon_computation_function_approx_(height);
	//     	rii_include::make_indices_convertion_function<double>(epsilon_grid, offset_f_)(output.data());    
 //    	}
 //    	else
 //    	{
 //    		auto epsilon_grid = epsilon_computation_function_(height);
	//     	rii_include::make_indices_convertion_function<double>(epsilon_grid, offset_f_)(output.data());    
 //    	}
    	
 //    	// set S_vec_ accordingly        	
 //    	ierr = VecSetValues(S_vec_, block_size, ix, &output[0], INSERT_VALUES);CHKERRV(ierr);		
 //    }
    
	// test S = I;
	// ierr = VecCopy(I_field, S_vec_);CHKERRV(ierr);		

	// test S = 1
	ierr = VecSet(RT_problem_->S_vec_, 1.0);CHKERRV(ierr);		

	// test S = 0;
	// ierr = VecZeroEntries(S_vec_);CHKERRV(ierr);		

	// // LTE S = W_T;
	// const auto block_size = RT_problem_->block_size_;
 		
	// int istart, iend; 
	// ierr = VecGetOwnershipRange(RT_problem_->I_field_, &istart, &iend);CHKERRV(ierr);	

	// size_t i_space;
	// double value;

	// for (int i = istart; i < iend; ++i)
	// {
	// 	i_space = i / block_size;

	// 	value = RT_problem_->W_T_[i_space];

	// 	ierr = VecSetValues(S_vec_, 1, &i, &value, INSERT_VALUES);CHKERRV(ierr); 							
	// }

	// ierr = VecAssemblyBegin(S_vec_);CHKERRV(ierr); 
 	// ierr = VecAssemblyEnd(S_vec_);CHKERRV(ierr); 
    
 // 	// scale emission by eta_I ----------------> TODO
	// int ixx[4];
	// double eta_I, eta_I_inv;	
	// double S[4];
		
	// for (int i = istart; i < iend; i = i + 4)
	// {
	// 	ierr = VecGetValues(RT_problem_->eta_field_, 1, &i, &eta_I);CHKERRV(ierr);	

	// 	eta_I_inv = 1.0 / eta_I;

	// 	std::iota(ixx, ixx + 4, i);

	// 	ierr = VecGetValues(S_vec_, 4, ixx, S);CHKERRV(ierr);	

	// 	S[0] = eta_I_inv * S[0];
	// 	S[1] = eta_I_inv * S[1];
	// 	S[2] = eta_I_inv * S[2];
	// 	S[3] = eta_I_inv * S[3];

	// 	ierr = VecSetValues(S_vec_, 4, ixx, S, INSERT_VALUES);CHKERRV(ierr); 									
	// }

	// ierr = VecAssemblyBegin(S_vec_);CHKERRV(ierr); 
 // 	ierr = VecAssemblyEnd(S_vec_);CHKERRV(ierr);  
}


void RT_solver::assemble_rhs(){ // TODO

  	if (mpi_rank_ == 0) std::cout << "\nAssembling right hand side..."; 

	PetscErrorCode ierr;

	const auto N_nu       = RT_problem_->N_nu_;
	const auto tot_size   = RT_problem_->tot_size_;
	const auto block_size = RT_problem_->block_size_;	
	const auto local_size = RT_problem_->local_size_;

	// allocate rhs vector 
	ierr = VecCreate(PETSC_COMM_WORLD, &rhs_);CHKERRV(ierr);    
	ierr = VecSetSizes(rhs_,local_size,tot_size);CHKERRV(ierr);   	
	ierr = VecSetFromOptions(rhs_);CHKERRV(ierr);

	ierr = VecSet(rhs_, 1.0);CHKERRV(ierr);		

	// // fill it
	// const auto N_nu       = RT_problem_->N_nu_;
	// // const auto eta_field  = RT_problem_->eta_field_;

	// Vec eps_th; // TODO use directly rhs_
	// double eta_i, eta_I, value;
	// size_t i_space;
	// int index_I, index_s_nu;

	// ierr = VecCreate(PETSC_COMM_WORLD, &eps_th);CHKERRV(ierr);    
	// ierr = VecSetSizes(eps_th,local_size,tot_size);CHKERRV(ierr);   	
	// ierr = VecSetFromOptions(eps_th);CHKERRV(ierr);

	// int istart, iend;
	// ierr = VecGetOwnershipRange(eps_th, &istart, &iend);CHKERRV(ierr);

	// // fill eps_th =  eps_c_th +  eps_l_th
	// for (int i = istart; i < iend; ++i)
	// {	
	// 	value = 0.0;

	// 	std::vector<size_t> local_idx = RT_problem_->local_to_global(i);
	// 	i_space = local_idx[0];

	// 	if (LTE_) // eps = W_T
	// 	{
	// 		if (local_idx[4] == 0) value = RT_problem_->W_T_[i_space]; // TODO divide by eta_I?
	// 	}
	// 	else // eps = (eps_c_th + eps_l_th_) / eta_I
	// 	{
	// 		// get eta_i		
	// 		ierr = VecGetValues(eta_field, 1, &i, &eta_i);
			
	// 		// first Stokes parameter
	// 		if (local_idx[4] == 0)				
	// 		{	
	// 			index_s_nu = N_nu * i_space + local_idx[3];
				
	// 			if (RT_problem_->enable_continuum_) 
	// 			{
	// 				// eps_c_th
	// 				value = (RT_problem_->eps_c_th_[index_s_nu]);	

	// 				// eps_l_th		
	// 				value += (RT_problem_->epsilon_[i_space]) * (RT_problem_->W_T_[i_space]) * (eta_i - RT_problem_->k_c_[index_s_nu]);				
	// 			}
	// 			else
	// 			{
	// 				// eps_l_th
	// 				value += (RT_problem_->epsilon_[i_space]) * (RT_problem_->W_T_[i_space]) * eta_i;				
	// 			}

	// 			// (eps_c_th + eps_l_th_) / eta_I
	// 			value /= eta_i;
	// 		} 
	// 		else
	// 		{
	// 			// get eta_I (!= eta_i)
	// 			index_I = i - local_idx[4];
	// 			ierr = VecGetValues(eta_field, 1, &index_I, &eta_I);

	// 			// eps_l_th / eta_i_l
	// 			value = eta_i * (RT_problem_->epsilon_[i_space]) * (RT_problem_->W_T_[i_space]) / eta_I;				
	// 		}			
	// 	}
		
	// 	ierr = VecSetValue(eps_th,i,value,INSERT_VALUES);CHKERRV(ierr);   
	// }

	// ierr = VecAssemblyBegin(eps_th);CHKERRV(ierr); 
	// ierr = VecAssemblyEnd(eps_th);CHKERRV(ierr);

	// // save_vec(eps_th, "../output/epsth.m" ,"eps_th"); 
		
	// // fill rhs_ from formal solve with bc
	// mf_ctx_.formal_solve(rhs_, eps_th, 1.0); 	
	
	// // clean
	// ierr = VecDestroy(&eps_th);CHKERRV(ierr);

	if (mpi_rank_ == 0) std::cout << "done" << std::endl;	
}


// matrix-free matrix vector multiplication y = (Id - LJ)x
PetscErrorCode UserMult(Mat mat,Vec x,Vec y){

	PetscErrorCode ierr; 

	void *ptr;
   	MatShellGetContext(mat,&ptr);
  	MF_context *mf_ctx_ = (MF_context *)ptr;

  	auto RT_problem = mf_ctx_->RT_problem_;

  	// compute new emission in S_field_ (or first S_vec_?)
  	mf_ctx_->update_emission(x);   

  	// copy to field type 
  	mf_ctx_->vec_to_field(RT_problem->S_field_, RT_problem->S_vec_); 

  	// fill rhs_ from formal solve with zero bc  	
	mf_ctx_->formal_solve_global(RT_problem->I_field_, RT_problem->S_field_, 0.0);	

	// copy intensity into petscvec format
	mf_ctx_->field_to_vec(RT_problem->I_field_, y);
  	
	// update I_out = I_in - I_fs (y = x - y)
	ierr = VecAYPX(y, -1.0, x);CHKERRQ(ierr);

  	return ierr;
}


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
