#include "RT_solver.hpp"
#include "sgrid_SliceHalo.hpp"


void MF_context::field_to_vec(const Field_ptr_t field, Vec &v)
{
	if (mpi_rank_ == 0) std::cout << "\nCopying field to Vec...";

	PetscErrorCode ierr; 
		
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

    if (mpi_rank_ == 0) std::cout << "done" << std::endl;
}


void MF_context::vec_to_field(Field_ptr_t field, const Vec &v)
{
	if (mpi_rank_ == 0) std::cout << "\nCopying Vec to field..." << std::endl;

	PetscErrorCode ierr; 
		
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

	if (mpi_rank_ == 0) std::cout << "\nApplying BC...";

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

    if (mpi_rank_ == 0) std::cout << "done" << std::endl;										
}


void MF_context::find_intersection(double theta, double chi, const double Z_down, const double Z_top, const double L, t_intersect *T) 
{
    // check theta
    if (theta == 0 or chi == 0) std::cout << "WARNING: ray direction not supported in find_intersection!" << std::endl;       
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


std::vector<t_intersect> MF_context::find_prolongation(double theta, double chi, const double dz, const double L) {

    //// some modifications 

    // check theta
    if (theta == 0 or chi == 0 or theta == PI/2 or chi == PI/2) std::cout << "WARNING: ray direction not supported!" << std::endl;       
    {        
        theta += 1e-16;
        chi   += 1e-16;                 
    }
    
    // check widths
    if (dz <= 0 ) std::cout << "WARNING: dz not positive" << std::endl;               
    if (L  <= 0 ) std::cout << "WARNING: L not positive"  << std::endl;       
    
    // unit vector in the direction of the ray (minus for different convection in formal solver)
    const double st = sin(theta);
    const double x = - st * cos(chi);
    const double y = - st * sin(chi); 
    const double z = - cos(theta);

    //////

    static int nmax = -1;
    
    std::vector<t_intersect> is;
    std::vector<t_xyinters> xyi;

    int dix = (x > 0. ? 1 : -1);
    int diy = (y > 0. ? 1 : -1);
    int diz = (z > 0. ? 1 : -1);

    double x1 = x/sqrt(x*x+y*y), y1 = y/sqrt(x*x+y*y);
    double z_inters = diz > 0 ? dz : -dz;

    double T = fabs(dz/z);

    // (1) if necessary, allocate sufficient memory for the intersections:
    int nmax1 = MAX(fabs(L/x1), fabs(L/y1)) + 1;

    if (nmax1 > nmax) 
    {
        nmax = nmax1;

        // is.reserve(nmax);
        // xyi.reserve(nmax);
    }

    // (2) find all the intersections
    double t = 0;
    double xg = 0, yg = 0;
    int i = 0;
    int ix=(dix>0 ? 0 : -1), iy=(diy>0 ? 0 : -1);
    do {
        double tx, ty, txy;

        t_xyinters xyi_tmp;

        if (dix<0) {
            tx = (L*ix-xg)/x1;
        }
        else {
            tx = (L*(ix+dix)-xg)/x1;
        }
        if (diy<0) {
            ty = (L*iy-yg)/y1;
        }
        else {
            ty = (L*(iy+diy)-yg)/y1;
        }
        
        if (tx<=ty) {
            txy = tx;
            xyi_tmp.plane = 0;
            ix += dix;
        }
        else {
            txy = ty;
            xyi_tmp.plane = 1;
            iy += diy;
        }
        xg += txy * x1;
        yg += txy * y1;
        t = fabs(sqrt(xg*xg+yg*yg)/st);
        if (t > T) { //last point in the xy plane?
            xg = xyi_tmp.x = x*T;
            yg = xyi_tmp.y = y*T;
            xyi_tmp.z = dz*diz;
            xyi_tmp.plane = 2;
            t = T;
        }
        else {
            xyi_tmp.x = xg;
            xyi_tmp.y = yg;
            xyi_tmp.z = t*z;
            if (xg<L*ix) ix--;
            if (yg<L*iy) iy--;
        }

        xyi.push_back(xyi_tmp);

        i++;
    } while(t != T);

    const int N = i;  

    // (3) set ix, iy, iz indices
    xyi[0].ix = xyi[0].iy = xyi[0].iz = 0;
    for (i=0; i<N; i++) {
        switch (xyi[i].plane) {
            case 0:
                xyi[i].ix = (i>0 ? xyi[i-1].ix+dix : dix);
                xyi[i].iy = (i>0 ? xyi[i-1].iy : 0);
                xyi[i].iz = (i>0 ? xyi[i-1].iz : 0);
                break;
            case 1:
                xyi[i].ix = (i>0 ? xyi[i-1].ix : 0);
                xyi[i].iy = (i>0 ? xyi[i-1].iy+diy : diy);
                xyi[i].iz = (i>0 ? xyi[i-1].iz : 0);
                break;
            case 2:
                xyi[i].ix = (i>0 ? xyi[i-1].ix : 0);
                xyi[i].iy = (i>0 ? xyi[i-1].iy : 0);
                xyi[i].iz = (i>0 ? xyi[i-1].iz+diz : diz);
                break;
        }
    }

    // (4) first N-1 points
    double u1=0, u2=0, v1=0, v2=0, u=0, v=0, norm;
    double x_last=0, y_last=0, z_last=0;

    t_intersect is_tmp;

    for (i=0; i<N-1; i++) {
        is_tmp.distance = sqrt((xyi[i].x-x_last)*(xyi[i].x-x_last) + (xyi[i].y-y_last)*(xyi[i].y-y_last) + (xyi[i].z-z_last)*(xyi[i].z-z_last));
        switch (xyi[i].plane) {
            case 0:
                for (int j=0; j<4; j++) is_tmp.ix[j] = xyi[i].ix;
                is_tmp.iy[0] = MIN(xyi[i].iy, xyi[i].iy+diy); is_tmp.iz[0] = MIN(0,diz);
                is_tmp.iy[1] = MAX(xyi[i].iy, xyi[i].iy+diy); is_tmp.iz[1] = is_tmp.iz[0];
                is_tmp.iy[2] = is_tmp.iy[1]; is_tmp.iz[2] = 0;  is_tmp.iz[2] = MAX(0,diz);
                is_tmp.iy[3] = is_tmp.iy[0]; is_tmp.iz[3] = 0;  is_tmp.iz[3] = is_tmp.iz[2];
                u1 = L * is_tmp.iy[0]; u2 = L * is_tmp.iy[1];
                v1 = MIN(0., z_inters); v2 = MAX(0., z_inters);
                u = xyi[i].y; v = xyi[i].z;
                break;
            case 1:
                for (int j=0; j<4; j++) is_tmp.iy[j] = xyi[i].iy;
                is_tmp.ix[0] = MIN(xyi[i].ix, xyi[i].ix+dix); is_tmp.iz[0] = MIN(0,diz);
                is_tmp.ix[1] = MAX(xyi[i].ix, xyi[i].ix+dix); is_tmp.iz[1] = is_tmp.iz[0];
                is_tmp.ix[2] = is_tmp.ix[1]; is_tmp.iz[2] = 0;  is_tmp.iz[2] = MAX(0,diz);
                is_tmp.ix[3] = is_tmp.ix[0]; is_tmp.iz[3] = 0;  is_tmp.iz[3] = is_tmp.iz[2];
                u1 = L * is_tmp.ix[0]; u2 = L * is_tmp.ix[1];
                v1 = MIN(0., z_inters); v2 = MAX(0., z_inters);
                u = xyi[i].x; v = xyi[i].z;
                break;
        }
        norm = 1.0 / ((u2-u1) * (v2-v1));
        is_tmp.w[0] = norm * (u2-u) * (v2-v);
        is_tmp.w[1] = norm * (u-u1) * (v2-v);
        is_tmp.w[2] = norm * (u-u1) * (v-v1);
        is_tmp.w[3] = norm * (u2-u) * (v-v1);
        x_last=xyi[i].x;
        y_last=xyi[i].y;
        z_last=xyi[i].z;

        is.push_back(is_tmp);
    }
    // the last point in the XY plane
    is_tmp.distance = sqrt((xyi[i].x-x_last)*(xyi[i].x-x_last) + (xyi[i].y-y_last)*(xyi[i].y-y_last) + (xyi[i].z-z_last)*(xyi[i].z-z_last));
    for (int j=0; j<4; j++) is_tmp.iz[j] = diz;
    is_tmp.ix[0] = MIN(xyi[i].ix, xyi[i].ix+dix); is_tmp.iy[0] = MIN(xyi[i].iy, xyi[i].iy+diy);
    is_tmp.ix[1] = MAX(xyi[i].ix, xyi[i].ix+dix); is_tmp.iy[1] = MIN(xyi[i].iy, xyi[i].iy+diy);
    is_tmp.ix[2] = MAX(xyi[i].ix, xyi[i].ix+dix); is_tmp.iy[2] = MAX(xyi[i].iy, xyi[i].iy+diy);
    is_tmp.ix[3] = MIN(xyi[i].ix, xyi[i].ix+dix); is_tmp.iy[3] = MAX(xyi[i].iy, xyi[i].iy+diy);
    u1 = L * is_tmp.ix[0]; u2 = L * is_tmp.ix[1];
    v1 = L * is_tmp.iy[0]; v2 = L * is_tmp.iy[3];
    u = xyi[i].x; v = xyi[i].y;
    norm = 1.0 / ((u2-u1) * (v2-v1));
    is_tmp.w[0] = norm * (u2-u) * (v2-v);
    is_tmp.w[1] = norm * (u-u1) * (v2-v);
    is_tmp.w[2] = norm * (u-u1) * (v-v1);
    is_tmp.w[3] = norm * (u2-u) * (v-v1);

    is.push_back(is_tmp);

    return is;
}

void MF_context::get_2D_weigths(const double x, const double y, double *w)
{
	// get weigths in the unit square
	if (x >= 1 or x <= 0) std::cout << "Problem in x input "<< std::endl;
	if (y >= 1 or y <= 0) std::cout << "Problem in y input "<< std::endl;

	const double xy = x * y;

	w[0] = 1.0 - x - y + xy; // (1.0 - x) * (1.0 - y) 
	w[1] = x - xy; // (1.0 - y) * x
	w[2] = xy;
	w[3] = y - xy;  //(1.0 - x) * y; 
}

// given a intersection type with N cells and grid indeces ijk, get I1, S1, K1 i.e. quantities needed for the last step of formal solution
std::vector<double> MF_context::long_ray_steps(const std::vector<t_intersect> T, 
                                               const Field_ptr_t I_field, const Field_ptr_t S_field, 
                                               const int i, const int j, const int k, const int block_index)
{     
    const auto N_x = RT_problem_->N_x_;
    const auto N_y = RT_problem_->N_y_;

    const auto eta_dev = RT_problem_->eta_field_->view_device();
    const auto rho_dev = RT_problem_->rho_field_->view_device();

    const auto I_dev = I_field->view_device();     
    const auto S_dev = S_field->view_device(); 

    const auto vertical_decomposition = RT_problem_->only_vertical_decomposition_;
    
	// coeff trap + cm conversion = - 0.5 * 1e5;
	const double coeff = -50000;

    // number of traversec cells 
    const int N = T.size();

	int i_intersect, j_intersect, k_intersect, b_index;

	double eta_I_1, weight, dtau;

	std::vector<double> I1(4), I2(4), S1(4), S2(4), etas(4), rhos(4), K1(16), K2(16);
   
	for (int cell = 0; cell < N; ++cell)
	{							
		// quantities in (1)
		if (cell == 0) 
		{
			// init
			for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
			{
				etas[i_stokes] = 0;
				rhos[i_stokes] = 0;
				S1[i_stokes]   = 0;
				I1[i_stokes]   = 0;
			}

			for (int face_vertices = 0; face_vertices < 4; ++face_vertices)
			{
				i_intersect = i + T[cell].ix[face_vertices];
				j_intersect = j + T[cell].iy[face_vertices];
				k_intersect = k - T[cell].iz[face_vertices]; 

				// correct for periodic boundary
				if (vertical_decomposition)
				{
					i_intersect = i_intersect % N_x; 
					j_intersect = j_intersect % N_y;
					if (i_intersect == 0) i_intersect = N_x;
					if (j_intersect == 0) j_intersect = N_y;										
				}	 																		

				weight = T[cell].w[face_vertices];	

				for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
				{
					b_index = block_index + i_stokes;										
			
					// get eta and rho
					etas[i_stokes] += weight * eta_dev.block(i_intersect,j_intersect,k_intersect)[b_index]; 
					rhos[i_stokes] += weight * rho_dev.block(i_intersect,j_intersect,k_intersect)[b_index];

					// set S1 and I1 
					S1[i_stokes] += weight * S_dev.block(i_intersect,j_intersect,k_intersect)[b_index];												
					I1[i_stokes] += weight * I_dev.block(i_intersect,j_intersect,k_intersect)[b_index];																									
				}	
			}

			K1 = assemble_propagation_matrix_scaled(etas, rhos);
		}
		else // reuse quantities in (2) 
		{
			S1 = S2;
			I1 = I2;
			K1 = K2;
		}        

		// save for later use
		eta_I_1 = etas[0];
		
		// quantities in (2)   
        if (cell < N - 1)
        {           
            // init
            for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
            {
                etas[i_stokes] = 0;
                rhos[i_stokes] = 0;
                S2[i_stokes]   = 0;            
            }

    		for (int face_vertices = 0; face_vertices < 4; ++face_vertices)
    		{
    			i_intersect = i + T[cell + 1].ix[face_vertices];
    			j_intersect = j + T[cell + 1].iy[face_vertices];
    			k_intersect = k - T[cell + 1].iz[face_vertices]; 

    			// correct for periodic boundary
    			if (vertical_decomposition)
    			{
    				i_intersect = i_intersect % N_x; 
    				j_intersect = j_intersect % N_y;
    				if (i_intersect == 0) i_intersect = N_x;
    				if (j_intersect == 0) j_intersect = N_y;										
    			}	 																		

    			weight = T[cell + 1].w[face_vertices];	

    			for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
    			{
    				b_index = block_index + i_stokes;										

    				// get eta and rho
    				etas[i_stokes] += weight * eta_dev.block(i_intersect,j_intersect,k_intersect)[b_index]; 
    				rhos[i_stokes] += weight * rho_dev.block(i_intersect,j_intersect,k_intersect)[b_index];

    				// set S2
    				S2[i_stokes] += weight * S_dev.block(i_intersect,j_intersect,k_intersect)[b_index];															
    			}				
    		}

    		K2 = assemble_propagation_matrix_scaled(etas, rhos);		
        }
        else
        {
            for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
            {
                b_index = block_index + i_stokes;                                       
        
                // get eta and rho
                etas[i_stokes] = eta_dev.block(i,j,k)[b_index]; 
                rhos[i_stokes] = rho_dev.block(i,j,k)[b_index];

                // set S2
                S2[i_stokes] = S_dev.block(i,j,k)[b_index];                                                         
            }               

            K2 = assemble_propagation_matrix_scaled(etas, rhos);
        }

		// optical depth step								
		dtau = coeff * (eta_I_1 + etas[0]) * T[cell].distance;									

		if (dtau > 0 ) std::cout << "ERROR in dtau sign" << std::endl;		
        
		formal_solver_.one_step(dtau, K1, K2, S1, S2, I1, I2);	
	}	
                                                                                                                    
    return I2;
}

void MF_context::formal_solve_local(Field_ptr_t I_field, const Field_ptr_t S_field, const Real I0)
{

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

	if (i_start > 1 or j_start > 1 or k_start > 1) std::cout << "WARNING: tested only for margin = 1!"<< std::endl;

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

	const auto I_dev = 	I_field->view_device();		
	const auto S_dev = 	S_field->view_device();	
	const auto g_dev = space_grid->view_device();

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	if (i_start > 1 or j_start > 1 or k_start > 1) std::cout << "WARNING: tested only for margin = 1!"<< std::endl;

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
	const int k_end = k_start + g_dev.dim[2];	

	int i_aux, j_aux, k_aux, k_global, i_intersect, j_intersect, k_intersect, b_start, b_index;

	// misc coeffs
	double theta, chi, mu, dtau, weight, eta_I_1, dz, cos_chi, sin_chi;	    
	
	bool boundary, horizontal_face, long_ray;

	// quantities depending on spatial point i
	std::vector<double> I1(4), I2(4), S1(4), S2(4), etas(4), rhos(4), K1(16), K2(16);

	// intersection object
	t_intersect intersection_data;        
 std::vector<t_intersect> intersection_data_long_ray;

	// minus for optical depth conversion, trap rule and conversion to cm (- 0.5 * 1e5)
	const double coeff = -50000;

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
							dz = (mu > 0) ? depth_grid[k_global] -  depth_grid[k_global + 1] : depth_grid[k_global - 1] - depth_grid[k_global];							

							for (int k_chi = 0; k_chi < (int)N_chi; ++k_chi)
							{	                                
								chi = chi_grid[k_chi];

								cos_chi = cos(chi);
								sin_chi = sin(chi);
							
								i_aux = (cos_chi < 0.0) ? i_end - i : i;	
								j_aux = (sin_chi < 0.0) ? j_end - j : j;		

								find_intersection(theta, chi, dz, dz, L, &intersection_data); 

								horizontal_face = intersection_data.iz[0] == intersection_data.iz[1] and 
									     		  intersection_data.iz[0] == intersection_data.iz[2] and 
									     		  intersection_data.iz[0] == intersection_data.iz[3];

                                // if long ray start in a different processor use short ray
								long_ray = (not horizontal_face) and (i == i_start or j == j_start) and (vertical_decomposition);								     		  	                                
                                								
								// check if a vertical face is intersected and use long ray
								if (long_ray)																																															
								{				                                        
                                    // intersection_data_long_ray = find_prolongation(theta, chi, dz, L);                                     
                                    intersection_data_long_ray = find_prolongation(theta, chi, dz, L);  
                                    
                                    // check number of cells traversed by the long ray
									if (intersection_data_long_ray.size() < 2) std::cout << "WARNING: number of traversed cells < 2 for long ray!"<< std::endl;
								}	


                                // set intersecion indeces
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
                                }

								// loop on freqs
								for (int n = 0; n < (int)N_nu; ++n) 
								{			                                    
									// block index
									b_start = RT_problem_->local_to_block(j_theta, k_chi, n); 

									// solve ODE
									if (long_ray)
									{                                    										
                                        I2 = long_ray_steps(intersection_data_long_ray, I_field, S_field, i_aux, j_aux, k_aux, b_start);
									}
									else // short ray
									{								
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

										// compute K1, S1 and I1						
																										
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
									}
									
									// test
									if (j_theta == N_theta/2 and k_chi == 0 and n == 0 and rank == mpi_size_ - 1)
									{									
										// std::cout << "I1 = " << I1[0] << std::endl;		
										if (i == 1 and j == 1 and k == k_end - 1)
										{
											// std::cout << intersection_data.distance << std::endl;	
											// std::cout << "coeff = " << coeff << std::endl;	
											// std::cout << "eta_I_1 = " << eta_I_1 << std::endl;	
											// std::cout << "etas[0] = " << etas[0] << std::endl;	
											// std::cout << "intersection_data.distance = " << intersection_data.distance << std::endl;					
											
											std::cout << "I1 = " << I1[0] << std::endl;	
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
		if (mpi_size_ > 1) I_field->exchange_halos(); // --------> TODO	better
	}
}

	
void MF_context::set_up_emission_module(){

    if (mpi_rank_ == 0) std::cout << "\nSetting up emission module...";
    
    // set some aliases 
    using rii_eps_comp_3D                 = rii_include::emission_coefficient_computation_3D;
    using rii_formal_solver_factory       = rii_include::formal_solver_factory_from_3D_RT_problem;
    using in_RT_problem_3D                = rii_include::in_RT_problem_3D_interface<RT_problem>;
    using emission_coefficient_components = rii_include::emission_coefficient_computation::emission_coefficient_components;

    // Build module
    ecc_sh_ptr_ = rii_eps_comp_3D::make_emission_coefficient_computation_3D_shared_ptr();

    auto fsf_sh_ptr = rii_formal_solver_factory::make_formal_solver_factory_from_3D_RT_problem_shared_ptr();

    in_RT_problem_3D::add_models(RT_problem_, ecc_sh_ptr_, fsf_sh_ptr, true);

    fsf_sh_ptr->make_formal_solver();

    std::list<emission_coefficient_components> components{
        emission_coefficient_components::epsilon_R_II,
        emission_coefficient_components::epsilon_R_III_GL,
        emission_coefficient_components::epsilon_csc};

    epsilon_fun_ = ecc_sh_ptr_->make_computation_function(components);

    offset_fun_ = rii_include::make_default_offset_function(RT_problem_->N_theta_, RT_problem_->N_chi_, RT_problem_->N_nu_);
        	
 //    epsilon_computation_function_approx_ = ecc_sh_ptr_->make_computation_function(
	// {
	// 	// rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_R_II,
	//  	rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_R_III_CRD_limit,
	//  	// rii_include::emission_coefficient_computation::emission_coefficient_components::epsilon_csc
	// },  
	// rii_consts::rii_units::kilometer);		  

    if (mpi_rank_ == 0) std::cout << "done" << std::endl;
}


// emissivity from scattering (line + continuum)
void MF_context::update_emission(const Vec &I_field, const bool approx){ 	

	PetscErrorCode ierr; 
	       
    const auto g_dev   = RT_problem_->space_grid_->view_device();  
    const auto eta_dev = RT_problem_->eta_field_->view_device();
    const auto S_dev   = RT_problem_->S_field_->view_device(); 

    // field range indeces 
    const int i_start = g_dev.margin[0]; 
    const int j_start = g_dev.margin[1];
    const int k_start = g_dev.margin[2];

    const int i_end = i_start + g_dev.dim[0];
    const int j_end = j_start + g_dev.dim[1];
    const int k_end = k_start + g_dev.dim[2];   

    const auto block_size = RT_problem_->block_size_;   
	
    std::vector<double>  input(block_size);        
    std::vector<double> output(block_size); 

    int ix[block_size];

    int istart, iend; 
    ierr = VecGetOwnershipRange(I_field, &istart, &iend);CHKERRV(ierr);	
	
	const int istart_local = istart / block_size;
	const int iend_local   = iend   / block_size;

    // grid indeces
    int i, j, k;    

    int counter_i = 0;
    int counter_j = 0;
    int counter_k = 0;
    
    for (int i_vec = istart_local; i_vec < iend_local; ++i_vec)
    {
    	// set indeces
    	std::iota(ix, ix + block_size, i_vec * block_size);

        // compute grid indeces from Vec index i_vec
        i = i_start + counter_i;
        j = j_start + counter_j;
        k = k_start + counter_k;

        if (i >= i_end or j >= j_end or k >= k_end) std::cout << "ERROR with counters in update_emission()!" << std::endl;
                
    	// get I field 
    	ierr = VecGetValues(I_field, block_size, ix, &input[0]);CHKERRV(ierr);	

    	// set input field
        ecc_sh_ptr_->update_incoming_field(i, j, k, offset_fun_, input.data());
       	
    	if (approx)
    	{
    		// auto epsilon_grid = epsilon_computation_function_approx_(height);
	    	// rii_include::make_indices_convertion_function<double>(epsilon_grid, offset_f_)(output.data());    
    	}
    	else
    	{
    		const auto out_field = epsilon_fun_(i,j,k);	    	
            rii_include::make_indices_convertion_function<double>(out_field, offset_fun_)(output.data());  
    	}
    	    	
        // update S_field_ from output scaling by eta_I
        double eta_I_inv; 

        for (int b = 0; b < block_size; b = b + 4)
        {
            eta_I_inv = 1.0 / (eta_dev.block(i,j,k)[b]);
            
            S_dev.block(i,j,k)[b]     = eta_I_inv * output[b];                                                         
            S_dev.block(i,j,k)[b + 1] = eta_I_inv * output[b + 1];
            S_dev.block(i,j,k)[b + 2] = eta_I_inv * output[b + 2];
            S_dev.block(i,j,k)[b + 3] = eta_I_inv * output[b + 3];                                                       
        }

        // update counters
        counter_i++;

        if (counter_i == i_end - 1)
        {
            counter_i = 0;
            counter_j++;
        }

        if (counter_j == j_end - 1)
        {
            counter_j = 0;
            counter_k++;
        }
    }

    // TODOs tetst
	// test S = I;
	// test S = 1
	// test S = 0;
}


void RT_solver::assemble_rhs(){ 

  	if (mpi_rank_ == 0) std::cout << "\n++++++ Assembling right hand side...+++++++++";
 
	PetscErrorCode ierr;

	const auto space_grid = RT_problem_->space_grid_;	

	// get sizes
	const auto N_nu       = RT_problem_->N_nu_;
	const auto tot_size   = RT_problem_->tot_size_;
	const auto block_size = RT_problem_->block_size_;	
	const auto local_size = RT_problem_->local_size_;

	// get fields
	const auto eta_dev =      RT_problem_->eta_field_->view_device();	
	const auto eps_c_th_dev = RT_problem_->eps_c_th_ ->view_device();	
	const auto epsilon_dev  = RT_problem_->epsilon_  ->view_device();	
	const auto W_T_dev      = RT_problem_->W_T_      ->view_device();	
	const auto k_c_dev      = RT_problem_->k_c_      ->view_device();	

	// allocate rhs vector 
	ierr = VecCreate(PETSC_COMM_WORLD, &rhs_);CHKERRV(ierr);    
	ierr = VecSetSizes(rhs_,local_size,tot_size);CHKERRV(ierr);   	
	ierr = VecSetFromOptions(rhs_);CHKERRV(ierr);

	// create rhs field (temporary)
	auto rhs_field = std::make_shared<Field_t>("EPS_TH", space_grid, block_size, sgrid::STAR_STENCIL);
	rhs_field->allocate_on_device(); 	
	
	// create eps_th field (temporary)
	auto eps_th_field = std::make_shared<Field_t>("EPS_TH", space_grid, block_size, sgrid::STAR_STENCIL);
	eps_th_field->allocate_on_device(); 
	const auto eps_th_dev = eps_th_field->view_device();	 
	
	// fill it eps_th =  eps_c_th +  eps_l_th
	sgrid::parallel_for("ASSEMBLE RHS", space_grid->md_range(), SGRID_LAMBDA(int i, int j, int k) {
		
		double value;

		std::vector<size_t> local_idx;

		for (int b = 0; b < (int)block_size; b++) 
		{		
			local_idx = RT_problem_->block_to_local(b);

			double eta_i = eta_dev.block(i,j,k)[b];

			// first Stokes parameter
			if (local_idx[3] == 0)
			{				
				auto index_nu = local_idx[2];
				
				if (RT_problem_->enable_continuum_) 
				{
					// eps_c_th
					value = eps_c_th_dev.block(i,j,k)[index_nu];	

					// eps_l_th		
					value += epsilon_dev.ref(i,j,k) * W_T_dev.ref(i,j,k) * (eta_i - k_c_dev.block(i,j,k)[index_nu]);				
				}
				else
				{
					// eps_l_th
					value = epsilon_dev.ref(i,j,k) * W_T_dev.ref(i,j,k) * eta_i;				
				}

				// (eps_c_th + eps_l_th_) / eta_I
				value /= eta_i;			
			}
			else
			{
				// get eta_I (!= eta_i)
				double eta_I = eta_dev.block(i,j,k)[b - local_idx[3]];

				// eps_l_th / eta_i_l
				value = eta_i * epsilon_dev.ref(i,j,k) * W_T_dev.ref(i,j,k) / eta_I;	
			}	

			// finally se eps_th
			eps_th_dev.block(i,j,k)[b] = value;	            
		}	
	});

	// fill rhs_ from formal solve with bc
	mf_ctx_.formal_solve_global(rhs_field, eps_th_field, 1.0); 	
	mf_ctx_.field_to_vec(rhs_field, rhs_);

	if (mpi_rank_ == 0) std::cout << "\n+++++++++++++++++++++++++++++++++++++++++++++\n"; 	
}


// matrix-free matrix vector multiplication y = (Id - LJ)x
PetscErrorCode UserMult(Mat mat,Vec x,Vec y){

	PetscErrorCode ierr; 

	void *ptr;
   	MatShellGetContext(mat,&ptr);
  	MF_context *mf_ctx_ = (MF_context *)ptr;

  	auto RT_problem = mf_ctx_->RT_problem_;

  	// compute new emission in S_field_ 
  	mf_ctx_->update_emission(x);   
  	
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
