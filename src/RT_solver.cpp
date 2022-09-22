#include "RT_solver.hpp"

//////////////////////////////////////////////////////
// Jiri functions for find_prolongation
static bool wError(const t_intersect &inters) {
    double ws = 0;
    for (int i=0; i<4; i++) {
        if (inters.w[i] < 0. || inters.w[i]>1.) return true;
        ws += inters.w[i];
    }
    if (fabs(ws-1.)>1e-10) return true;
    return false;
}

// r is "almost an integer" with possible small rounding error; convert it to int
inline static int r2int(double r) {
    return r >= 0. ? (int)(r+1e-3) : (int)(r-1e-3);
}

static void setIXYZ(double L, t_xyzinters &intersect) {
    switch (intersect.plane) {
        case I_YZ:
            intersect.ix = r2int(intersect.x/L);
            intersect.iy = intersect.y>=0. ? int(intersect.y/L) : int(intersect.y/L)-1;
            intersect.iz = intersect.z>0. ? 0 : -1;
            break;
        case I_XZ:
            intersect.ix = intersect.x>=0. ? int(intersect.x/L) : int(intersect.x/L)-1;
            intersect.iy = r2int(intersect.y/L);
            intersect.iz = intersect.z>0. ? 0 : -1;
            break;
        case I_XY:
            intersect.ix = intersect.x>=0. ? int(intersect.x/L) : int(intersect.x/L)-1;
            intersect.iy = intersect.y>=0. ? int(intersect.y/L) : int(intersect.y/L)-1;
            intersect.iz = intersect.z>0. ? 1 : -1;
            break;
        default:
            std::cout << "WARNING: Invalid plane!" << std::endl;
    }
}


bool intersComp(const t_xyzinters &i1, const t_xyzinters &i2) {
    return (i1.t < i2.t);
}
//////////////////////////////////////////////////////

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

	int istart, row;

	double value;
	
	ierr = VecGetOwnershipRange(v, &istart, NULL);CHKERRV(ierr);	

	int counter = 0;

	for (int k = k_start; k < k_end; ++k)					
	{															
		for (int j = j_start; j < j_end; ++j)
		{
			for (int i = i_start; i < i_end; ++i)				
			{
				for (int b = 0; b < (int)block_size; b++) 
				{
					// set row index and correposnding entry
					row = istart + counter;

					value = f_dev.block(i, j, k)[b];

					ierr = VecSetValue(v, row, value, INSERT_VALUES);CHKERRV(ierr); // TODO: use VecSetValues for perf

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

	int istart, row;

	double value;
	
	ierr = VecGetOwnershipRange(v, &istart, NULL);CHKERRV(ierr);	

	int counter = 0;

	for (int k = k_start; k < k_end; ++k)					
	{															
		for (int j = j_start; j < j_end; ++j)
		{
			for (int i = i_start; i < i_end; ++i)				
			{
				for (int b = 0; b < (int)block_size; b++) 
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

    // init some quantities 
    const auto N_z        = RT_problem_->N_z_;
    const auto block_size = RT_problem_->block_size_;
    const auto space_grid = RT_problem_->space_grid_;   

    // apply BC
    const auto W_T_dev     = RT_problem_->W_T_->view_device();
    const auto g_dev       =        space_grid->view_device();
    auto I_field_dev       =           I_field->view_device();        

    sgrid::parallel_for("APPLY BC", space_grid->md_range(), SGRID_LAMBDA(int i, int j, int k) {
                                    
        // just in max depth
        if (g_dev.global_coord(2, k) == (N_z - 1))         
        {       
            const Real W_T_deep = I0 * W_T_dev.ref(i,j,k);
            
            for (int b = 0; b < (int)block_size; b = b + 4) 
            {
                I_field_dev.block(i,j,k)[b] = W_T_deep;         
            }            
        }
    }); 
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

    // sanity checks
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

        const double w_sum = T->w[0] + T->w[1] + T->w[2] + T->w[3];

        if (std::abs(w_sum - 1.0) > 1e-15) std::cout << "WARNING: w_sum - 1 = " << w_sum - 1.0 << ", it should be 0!" << std::endl;            
    }
}

////////////////////////////////////
#define DEBUG_MODE
////////////////////////////////////


std::vector<t_intersect> MF_context::find_prolongation(double theta, double chi, const double dz, const double L) {

    #ifdef DEBUG_MODE
    const double small = 1e-10;
    if (theta<0. || theta>PI || chi<0. || chi>2.*PI) {
        std::cout << "WARNING: angles out of allowed intervals!" << std::endl;
    }
    if (fabs(theta)<small || fabs(theta-PI/2)<small || fabs(theta-PI)<small ||
        fabs(chi)<small || fabs(chi-PI/2)<small || fabs(chi-PI)<small || fabs(chi-3.*PI/2)<small || fabs(chi-2.*PI)<small) {
        std::cout << "WARNING: ray direction not supported!" << std::endl;
        theta += small;
        chi   += small;                 
    }

    // check widths
    if (dz <= 0 ) std::cout << "WARNING: dz not positive" << std::endl;
    if (L  <= 0 ) std::cout << "WARNING: L not positive"  << std::endl;
    #endif
    
    // unit vector in the direction of the ray (minus for different convention in formal solver)
    const double st = sin(theta);
    const double x = - st * cos(chi);
    const double y = - st * sin(chi); 
    const double z = - cos(theta);

    std::vector<t_intersect> is;
    std::vector<t_xyzinters> xyi;

    // step directions:
    int dix = (x > 0. ? 1 : -1);
    int diy = (y > 0. ? 1 : -1);
    int diz = (z > 0. ? 1 : -1);

    double x1 = x/sqrt(x*x+y*y), y1 = y/sqrt(x*x+y*y); // unit vector of the ray projection to the xy plane
    double T = fabs(dz/z); // total length of the ray

    // ALGORITHM: find all intersections with x and y planes, add them to a vector and sort w/ respect to distance:
    t_xyzinters xyi_tmp;

    // 1) add the final intersection with the xy plane
    xyi_tmp.plane = I_XY;
    xyi_tmp.x = x*T;    xyi_tmp.y = y*T;    xyi_tmp.z = z*T;
    xyi_tmp.t = T;
    setIXYZ(L, xyi_tmp);
    xyi.push_back(xyi_tmp); // add to the list as the first element

    // 2) add the intersections with the xz planes:
    double t = 0;
    for (int iy=diy; t<T; iy+=diy) {
        double txy = L*iy / y1;
        xyi_tmp.x = x1*txy;
        xyi_tmp.y = L*iy;
        xyi_tmp.t = t = fabs(sqrt(xyi_tmp.x*xyi_tmp.x+xyi_tmp.y*xyi_tmp.y)/st);
        xyi_tmp.z = z*t;
        xyi_tmp.plane = I_XZ;
        setIXYZ(L, xyi_tmp);
        if (t<T) {
            xyi.push_back(xyi_tmp);
        }
    }
    
    // 3) add the intersections with the yz planes:
    t = 0;
    for (int ix=dix; t<T; ix+=dix) {
        double txy = L*ix / x1;
        xyi_tmp.x = L*ix;
        xyi_tmp.y = y1*txy;
        xyi_tmp.t = t = fabs(sqrt(xyi_tmp.x*xyi_tmp.x+xyi_tmp.y*xyi_tmp.y)/st);
        xyi_tmp.z = z*t;
        xyi_tmp.plane = I_YZ;
        setIXYZ(L, xyi_tmp);
        if (t<T) {
            xyi.push_back(xyi_tmp);
        }
    }

    // 4) sort the vector w/ respect to t
    sort(xyi.begin(), xyi.end(), intersComp);

    #ifdef DEBUG_MODE
    if (xyi[xyi.size()-1].plane != I_XY) std::cout << "WARNING: the last intersation is not in the XY plane!" << std::endl;
    #endif

    for (unsigned int i=0; i<xyi.size(); i++) {
        t_intersect is_tmp;
        double u1, u2, v1, v2, u, v, norm;

        is_tmp.distance = xyi[i].t;
        switch (xyi[i].plane) {
            case I_YZ:
                for (int j=0; j<4; j++) is_tmp.ix[j] = xyi[i].ix;
                is_tmp.iy[0] = xyi[i].iy;   is_tmp.iz[0] = xyi[i].iz;
                is_tmp.iy[1] = xyi[i].iy+1; is_tmp.iz[1] = xyi[i].iz;
                is_tmp.iy[2] = xyi[i].iy+1; is_tmp.iz[2] = xyi[i].iz+1;
                is_tmp.iy[3] = xyi[i].iy;   is_tmp.iz[3] = xyi[i].iz+1;
                u1 = L * is_tmp.iy[0]; u2 = L * is_tmp.iy[1];
                v1 = diz > 0 ? 0. : -dz;  v2 = diz > 0 ? dz : 0.;
                u = xyi[i].y; v = xyi[i].z;
                break;
            case I_XZ:
                for (int j=0; j<4; j++) is_tmp.iy[j] = xyi[i].iy;
                is_tmp.ix[0] = xyi[i].ix;   is_tmp.iz[0] = xyi[i].iz;
                is_tmp.ix[1] = xyi[i].ix+1; is_tmp.iz[1] = xyi[i].iz;
                is_tmp.ix[2] = xyi[i].ix+1; is_tmp.iz[2] = xyi[i].iz+1;
                is_tmp.ix[3] = xyi[i].ix;   is_tmp.iz[3] = xyi[i].iz+1;
                u1 = L * is_tmp.ix[0];    u2 = L * is_tmp.ix[1];
                v1 = diz > 0 ? 0. : -dz;  v2 = diz > 0 ? dz : 0.;
                u = xyi[i].x; v = xyi[i].z;
                break;
            case I_XY:
                for (int j=0; j<4; j++) is_tmp.iz[j] = diz > 0 ? 1 : -1;
                is_tmp.ix[0] = xyi[i].ix;   is_tmp.iy[0] = xyi[i].iy;
                is_tmp.ix[1] = xyi[i].ix+1; is_tmp.iy[1] = xyi[i].iy;
                is_tmp.ix[2] = xyi[i].ix+1; is_tmp.iy[2] = xyi[i].iy+1;
                is_tmp.ix[3] = xyi[i].ix;   is_tmp.iy[3] = xyi[i].iy+1;
                u1 = L * is_tmp.ix[0]; u2 = L * is_tmp.ix[1];
                v1 = L * is_tmp.iy[0]; v2 = L * is_tmp.iy[3];
                u = xyi[i].x; v = xyi[i].y;
                break;
            default:
                u1 = v1 = u2 = v2 = u = v = 0;
                std::cout << "WARNING: inconsistent intersection!" << std::endl;
        }
        norm = 1.0 / ((u2-u1) * (v2-v1));
        is_tmp.w[0] = norm * (u2-u) * (v2-v);
        is_tmp.w[1] = norm * (u-u1) * (v2-v);
        is_tmp.w[2] = norm * (u-u1) * (v-v1);
        is_tmp.w[3] = norm * (u2-u) * (v-v1);

        #ifdef DEBUG_MODE
        // sanity check
        if (wError(is_tmp)) {
            std::cout << "WARNING: w has a problem!" << std::endl;          
            std::cout << "theta = " << theta << std::endl;          
            std::cout << "chi = "  << chi << std::endl;             
            std::cout << "dz = "  << dz << std::endl;               
            std::cout << "L = "  << L << std::endl;             
                
            std::cout << "i = " << i << std::endl;          
            for (int k=0; k<4; k++) std::cout << "  w[" << k << "] = " << is_tmp.w[k] << std::endl;
        }
        #endif

        is.push_back(is_tmp);
    }

    // put initial condition in the first cell         
    std::reverse(is.begin(),is.end()); 
    
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

    std::cout << "WARNING: update distances in find prolongation! now buggy!" << std::endl; // TODO

    if (use_log_interpolation_) std::cout << "WARNING: log_interpolation not suppoerted for long_ray_steps()!" << std::endl;

    const auto N_x = RT_problem_->N_x_;
    const auto N_y = RT_problem_->N_y_;
    
    const auto eta_dev = eta_field_serial_->view_device(); 
    const auto rho_dev = rho_field_serial_->view_device(); 

    const auto I_dev = I_field->view_device();     
    const auto S_dev = S_field->view_device(); 
    
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

			    // correction for periodic BC 
                i_intersect = i_intersect % N_x; 
                j_intersect = j_intersect % N_y;            

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

    			// correction for periodic boundary
                i_intersect = i_intersect % N_x;
                j_intersect = j_intersect % N_y;

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


// given a intersection type with N cells and grid indeces ijk, get I1, S1, K1 i.e. quantities needed for the last step of formal solution
std::vector<double> MF_context::single_long_ray_step(const std::vector<t_intersect> T, 
                                               const Field_ptr_t I_field, const Field_ptr_t S_field, 
                                               const int i, const int j, const int k, const int block_index)
{             
    const auto N_x = RT_problem_->N_x_;
    const auto N_y = RT_problem_->N_y_;

    const auto eta_dev = eta_field_serial_->view_device(); 
    const auto rho_dev = rho_field_serial_->view_device(); 

    const auto I_dev = I_field->view_device();     
    const auto S_dev = S_field->view_device(); 

    // coeff trap + cm conversion = - 0.5 * 1e5;
    const double coeff = -50000;
    
    int i_intersect, j_intersect, k_intersect, b_index;

    double eta_I_1, weight;
    double total_distance = 0;

    std::vector<double> I1(4), I2(4), S1(4), S2(4), etas(4), rhos(4), K1(16), K2(16);
   
    // // compute total distance 
    // for (int cell = 0; cell < T.size(); ++cell)
    // {                           
    //    total_distance += T[cell].distance;
    // }   

    // TODO
    total_distance = T[0].distance;   // ----------------> FIXME in new version

    // quantities in (1)  
    for (int i_stokes = 0; i_stokes < 4; ++i_stokes) 
    {
        if (use_log_interpolation_)
        {
            etas[i_stokes] = 1.0;
            rhos[i_stokes] = 1.0;
            S1[i_stokes]   = 1.0;
            I1[i_stokes]   = 1.0;                                           
        }
        else // linear
        {
            etas[i_stokes] = 0;
            rhos[i_stokes] = 0;
            S1[i_stokes]   = 0;
            I1[i_stokes]   = 0;
        }
    }

    const double debug_value = std::abs(T[0].iz[0] + T[0].iz[1] + T[0].iz[2] + T[0].iz[3]);

    if (debug_value != 4) std::cout << "ERROR in single_long_ray_step()" << std::endl;

    for (int face_vertices = 0; face_vertices < 4; ++face_vertices)
    {
        i_intersect = i + T[0].ix[face_vertices];
        j_intersect = j + T[0].iy[face_vertices];
        k_intersect = k - T[0].iz[face_vertices]; 
        
        // correction for periodic BC 
        i_intersect = i_intersect % N_x; 
        j_intersect = j_intersect % N_y;
       
        weight = T[0].w[face_vertices];  

        if (weight < 0) std::cout << "weight = " << weight << std::endl;      
        
        for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
        {
            b_index = block_index + i_stokes;                                       

            if (use_log_interpolation_)
            {
                // get eta and rho                                                                                                    
                etas[i_stokes] *= pow_gen(eta_dev.block(i_intersect,j_intersect,k_intersect)[b_index], weight);                                                     
                rhos[i_stokes] *= pow_gen(rho_dev.block(i_intersect,j_intersect,k_intersect)[b_index], weight);
                
                // set S1 and I1
                S1[i_stokes] *= pow_gen(S_dev.block(i_intersect,j_intersect,k_intersect)[b_index], weight);
                I1[i_stokes] *= pow_gen(I_dev.block(i_intersect,j_intersect,k_intersect)[b_index], weight);
                
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

    // save for later use
    eta_I_1 = etas[0];

    if (eta_I_1 < 0 ) std::cout << "etas[0] = " << etas[0] << std::endl;      

    // quantities in (2)
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

    // optical depth step                               
    const double dtau = coeff * (eta_I_1 + etas[0]) * total_distance;                                  

    if (dtau > 0 ) std::cout << "ERROR in dtau sign" << std::endl;          
    
    formal_solver_.one_step(dtau, K1, K2, S1, S2, I1, I2);    
                                                                                                                            
    return I2;
}


// TODO: Why do I need serial fields as input? Use only I_field_serial_??

void MF_context::formal_solve_global(Field_ptr_t I_field, Field_ptr_t I_field_serial, 
                                     const Field_ptr_t S_field, const Field_ptr_t S_field_serial, const Real I0)
{
	if (mpi_rank_ == 0) std::cout << "\nStart global formal solution..." << std::endl;
    
	// init some quantities 	    
    const auto N_x = RT_problem_->N_x_;
    const auto N_y = RT_problem_->N_y_;
    const auto N_z = RT_problem_->N_z_;

    const auto N_theta = RT_problem_->N_theta_;
	
	const auto block_size = RT_problem_->block_size_;
	const auto tot_size   = RT_problem_->tot_size_;
	
	const auto mu_grid    = RT_problem_->mu_grid_;
	const auto theta_grid = RT_problem_->theta_grid_;	
	const auto chi_grid   = RT_problem_->chi_grid_;	
	const auto depth_grid = RT_problem_->depth_grid_;	
	const auto L          = RT_problem_->L_;		

    const auto eta_dev = eta_field_serial_->view_device(); 
    const auto rho_dev = rho_field_serial_->view_device(); 

    const auto I_dev = I_field_serial->view_device();		
	const auto S_dev = S_field_serial->view_device();	    

	const auto g_dev = space_grid_serial_->view_device();   

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	if (i_start > 0 or j_start > 0 or k_start > 0) std::cout << "WARNING: periodic BC hardcoded for margin = 0!" << std::endl;

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
	const int k_end = k_start + g_dev.dim[2];	
    
	int i_aux, j_aux, k_aux, k_global, b_start, b_index;

    std::vector<int> i_intersect(4), j_intersect(4), k_intersect(4);

    // serial indexing coeffs    
    std::vector<size_t> local_idx;
    int block_start, block_end, j_theta_start, k_chi_start, n_nu_start, j_theta_end, k_chi_end, n_nu_end;

	// misc coeffs
	double theta, chi, mu, dtau, weight, eta_I_1, dz;
	
	bool boundary, horizontal_face, long_ray;

	// quantities depending on spatial point i
	std::vector<double> I1(4), I2(4), S1(4), S2(4), etas(4), rhos(4), K1(16), K2(16);

	// intersection object
	t_intersect intersection_data;        
    std::vector<t_intersect> intersection_data_long_ray;

	// minus for optical depth conversion, trap rule and conversion to cm (- 0.5 * 1e5)
	const double coeff = -50000;
	
    double comm_timer     = 0;
    double one_step_timer = 0;
    double total_timer    = 0;

    const bool timing_debug = false;

    if (timing_debug) MPI_Barrier(MPI_COMM_WORLD);
    Real start_total = MPI_Wtime();                                    
    
    // impose boundary conditions 
    apply_bc(I_field, I0);  
    	
    for (int tile_number = 0; tile_number < n_tiles_; ++tile_number)  
	{	
        // get local block range
        block_start = mpi_rank_ * n_local_rays_ + tile_number * tile_size_; 
        block_end   = block_start + tile_size_ - 1;
        
        const bool idle_proc = (block_start > block_size - 1);
                         
        // communication	                
        if (timing_debug) MPI_Barrier(MPI_COMM_WORLD);
        Real start_comm = MPI_Wtime();                                    
        
        // write S to the serial grid and I to get initial condition              // TODO: why S_field_serial_?
        S_remap_.from_pgrid_to_pblock(*S_field, *S_field_serial_, tile_number);                        
        I_remap_.from_pgrid_to_pblock(*I_field, *I_field_serial_, tile_number); // TODO: this is a bit redundant, only one xy plane is needed
                
        comm_timer += MPI_Wtime() - start_comm;      

        if (not idle_proc)
        {      
            // get indeces
            local_idx = RT_problem_->block_to_local(block_start);
        
            j_theta_start = local_idx[0];
            k_chi_start   = local_idx[1];
            n_nu_start    = local_idx[2];

            bool throw_error = false;

            if (local_idx[3] != 0) { std::cout << "ERROR in block decomposition in formal_solve_global(), i_stokes_start not 0" << std::endl; throw_error = true; }

            local_idx = RT_problem_->block_to_local(block_end);

            j_theta_end = local_idx[0] + 1;
            k_chi_end   = local_idx[1] + 1;
            n_nu_end    = local_idx[2] + 1;                

            if (j_theta_end < j_theta_start) { std::cout << "ERROR with j_theta partition" << std::endl; throw_error = true; }      
            if (k_chi_end   < k_chi_start)   { std::cout << "ERROR with k_chi partition"   << std::endl; throw_error = true; }     
            if (n_nu_end    < n_nu_start)    { std::cout << "ERROR with n_nu partition"    << std::endl; throw_error = true; }   
        
            if (local_idx[3] != 3) { std::cout << "ERROR in block decomposition in formal_solve_global(), i_stokes_end not 3" << std::endl; throw_error = true; }      
            
            if (throw_error) throw "ERROR with block decomposition";
        
    		// loop over spatial points
    		for (int k = k_start; k < k_end; ++k)					
    		{									            
    			for (int j = j_start; j < j_end; ++j)
    			{
    				for (int i = i_start; i < i_end; ++i)
    				{					                       
    					// loop over directions (TODO could be parallel)
    					for (int j_theta = j_theta_start; j_theta < (int)j_theta_end; ++j_theta)
    					{
    						theta = theta_grid[j_theta];
    						mu    = mu_grid[j_theta];						

    						k_aux = (mu > 0.0) ? k_end - 1 - k + g_dev.margin[2]: k; 

    						// depth index
    						k_global = g_dev.global_coord(2, k_aux);	

    						boundary = (k_global == 0 and mu < 0) or (k_global == (int)N_z - 1 and mu > 0);
    						
    						if (not boundary)
    						{						
    							// set vertical box size
    							dz = (mu > 0) ? depth_grid[k_global] -  depth_grid[k_global + 1] : depth_grid[k_global - 1] - depth_grid[k_global];							

    							for (int k_chi = k_chi_start; k_chi < (int)k_chi_end; ++k_chi)
    							{	                                
    								chi = chi_grid[k_chi];
    								
    								i_aux = (cos(chi) < 0.0) ? i_end - 1 - i + g_dev.margin[0]: i;	
    								j_aux = (sin(chi) < 0.0) ? j_end - 1 - j + g_dev.margin[1]: j;		                                
                                    
    								find_intersection(theta, chi, dz, dz, L, &intersection_data); 

    								horizontal_face = intersection_data.iz[0] == intersection_data.iz[1] and 
    									     		  intersection_data.iz[0] == intersection_data.iz[2] and 
    									     		  intersection_data.iz[0] == intersection_data.iz[3];

                                    // menage short/long ray                               
                                    if (use_always_long_ray_)
                                    {
                                        long_ray = not horizontal_face;
                                    }
                                    else
                                    {
                                        long_ray = (not horizontal_face) and (i == i_start or j == j_start);     
                                    }     
                                        
                                    if (long_ray) intersection_data_long_ray = find_prolongation(theta, chi, dz, L);  
                                    
                                    // set intersection indeces
                                    if (not long_ray)
                                    {
                                        for (int face_v = 0; face_v < 4; ++face_v)
                                        {
                                            i_intersect[face_v] = i_aux + intersection_data.ix[face_v];
                                            j_intersect[face_v] = j_aux + intersection_data.iy[face_v];
                                            k_intersect[face_v] = k_aux - intersection_data.iz[face_v]; // minus because k increases going downwards  
                                            
                                            // impose periodic BC
                                            i_intersect[face_v] = i_intersect[face_v] % N_x;
                                            j_intersect[face_v] = j_intersect[face_v] % N_y;                                        
                                        }                                                           
                                    }                                                                 
                                    
    								// loop on freqs
    								for (int n = n_nu_start; n < (int)n_nu_end; ++n) 
    								{			     
                                        if (timing_debug) MPI_Barrier(MPI_COMM_WORLD);                                                                                                                 
                                        Real start_one = MPI_Wtime();                                               

    									// block index (coorected for tile size)                                    
    									b_start = RT_problem_->local_to_block(j_theta, k_chi, n) % tile_size_;                                   

    									// solve ODE
    									if (long_ray)
    									{                                                     
                                            if (use_single_long_step_)
                                            {
                                                I2 = single_long_ray_step(intersection_data_long_ray, I_field_serial, S_field_serial, i_aux, j_aux, k_aux, b_start);
                                            }
                                            else
                                            {
                                                I2 = long_ray_steps(intersection_data_long_ray, I_field_serial, S_field_serial, i_aux, j_aux, k_aux, b_start);
                                            }                                                                         
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
    												etas[i_stokes] = 1.0;
    												rhos[i_stokes] = 1.0;
    												S1[i_stokes]   = 1.0;
    												I1[i_stokes]   = 1.0;											
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
    										for (int face_v = 0; face_v < 4; ++face_v)
    										{			                                            
                                                weight = intersection_data.w[face_v];
    										
    											for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
    											{
    												b_index = b_start + i_stokes;										

    												if (use_log_interpolation_)
    												{
                                                        // get eta and rho                                                                                                    
                                                        etas[i_stokes] *= pow_gen(eta_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index], weight);                                                     
                                                        rhos[i_stokes] *= pow_gen(rho_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index], weight);
                                                        
                                                        // set S1 and I1
                                                        S1[i_stokes] *= pow_gen(S_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index], weight);
                                                        I1[i_stokes] *= pow_gen(I_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index], weight);
                                                        
    												}
    												else
    												{
    													// get eta and rho
    													etas[i_stokes] += weight * eta_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index]; 
    													rhos[i_stokes] += weight * rho_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index];                                                    

    													// set S1 and I1
    													S1[i_stokes] += weight * S_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index];	                                                    
    													I1[i_stokes] += weight * I_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index];	                                                                                                           
    												}									
    											}											
    										}																											                                        

    										K1 = assemble_propagation_matrix_scaled(etas, rhos);                                                                            
    										
    										// optical depth step								
    										dtau = coeff * (eta_I_1 + etas[0]) * intersection_data.distance;									
                                            
    										if (dtau > 0 ) std::cout << "ERROR in dtau sign" << std::endl;										
                                            
    										formal_solver_.one_step(dtau, K1, K2, S1, S2, I1, I2);	                                                                            
    									}
                                        
                                        one_step_timer += MPI_Wtime() - start_one;                       
    									                                    
    									// // test
    									// if (j_theta >= N_theta/2 and k_chi == 0 and n == 98)
    									// {									                                        
             //                                if (long_ray) std::cout << "WARNING LONG RAY: look at long ray routines for data!" << std::endl;
        //                                
                                            
             //                                std::cout << "theta = " << theta << std::endl;
             //                                // std::cout << "chi = "<< chi << std::endl;
             //                                std::cout << "mu = " << mu << std::endl;
             //                                // std::cout << "n = "  << n << std::endl;
                                             
    									// 	// std::cout << "I1 = " << I1[0] << std::endl;		
    									// 	// if (k == k_end - 1)
             //                                // if (i == 1 and j == 1)
    									// 	{											
             //                                    // std::cout << "dz = "<< dz << std::endl;

    									// 		// std::cout << "coeff = " << coeff << std::endl;	
    									// 		// std::cout << "eta_I_1 = " << eta_I_1 << std::endl;	
    									// 		// std::cout << "etas[0] = " << etas[0] << std::endl;	
    									// 		// std::cout << "intersection_data.distance = " << intersection_data.distance << std::endl;																                                            

             //                                    // std::cout << "mpi_rank_ = " << mpi_rank_ << std::endl;   
             //                                    std::cout << "k_global = " << g_dev.global_coord(2, k_aux) << std::endl;                                              
             //                                    // std::cout << "k = " << k << std::endl;                                              

             // 									std::cout << "I1 = "   << I1[0] << std::endl;	
             //                                    // std::cout << "Q1 = "   << I1[1] << std::endl;   
             //                                    // std::cout << "U1 = "   << I1[2] << std::endl;   
             //                                    // std::cout << "V1 = "   << I1[3] << std::endl;   

    									//  		std::cout << "I2 = "   << I2[0] << std::endl;    
             //                                    // std::cout << "Q2 = "   << I2[1] << std::endl;   
             //                                    // std::cout << "U2 = "   << I2[2] << std::endl;   
             //                                    // std::cout << "V2 = "   << I2[3] << std::endl;   

             //                                    // std::cout << "dtau = " << dtau  << std::endl; 	
    									// 		// std::cout << "dz = " << dz << std::endl;					    													
    									// 		// std::cout << "L = "  << L  << std::endl;					    											
    									// 		// std::cout << "chi = "  << chi  << std::endl;	
    											
    									//  	// 	std::cout << "S1 = " << std::endl;
    									//  	// 	for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S1[i_stokes] << std::endl;

    									//  	// 	std::cout << "S2 = " << std::endl;
    									//  	// 	for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S2[i_stokes] << std::endl;
    																		
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

        }      
          
        if (timing_debug) MPI_Barrier(MPI_COMM_WORLD);
        start_comm = MPI_Wtime();    
        
        I_remap_.from_pblock_to_pgrid(*I_field_serial, *I_field, tile_number); 

        comm_timer += MPI_Wtime() - start_comm; 
    }

    total_timer += MPI_Wtime() - start_total;                     
    
    if (mpi_rank_ == 0)
    {
        printf("comm_timer:\t\t%g seconds\n", comm_timer);
        printf("one_step_timer:\t\t%g seconds\n", one_step_timer);                        
        printf("total_timer:\t\t%g seconds\n", total_timer);                        
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
        // emission_coefficient_components::epsilon_R_II,
        emission_coefficient_components::epsilon_R_II_CONTRIB, // less memory is used for this version
        // emission_coefficient_components::epsilon_R_III,
        emission_coefficient_components::epsilon_R_III_GL,
        emission_coefficient_components::epsilon_csc};

    std::list<emission_coefficient_components> components_approx{
        // emission_coefficient_components::epsilon_R_II,
        emission_coefficient_components::epsilon_R_III_pCRD_limit,
        emission_coefficient_components::epsilon_csc
    };
    
    epsilon_fun_        = ecc_sh_ptr_->make_computation_function(components);    
    // if (mpi_rank_ == 0) std::cout << "================ emission TEST: solo epsilon_R_III_CRD_limit! ================" << std::endl;
    // epsilon_fun_        = ecc_sh_ptr_->make_computation_function(components_approx); // test
    
    epsilon_fun_approx_ = ecc_sh_ptr_->make_computation_function(components_approx);

    offset_fun_ = rii_include::make_default_offset_function(RT_problem_->N_theta_, RT_problem_->N_chi_, RT_problem_->N_nu_);
        	
	// rii_consts::rii_units::kilometer);		  //

    // set threads number
    // ecc_sh_ptr_->set_threads_number(2);

    // print memory
    unsigned long long b;
    b = ecc_sh_ptr_->bytes();

    if (mpi_rank_ == 0) std::cout << "\n[Memory from set_up_emission_module() = " << (double)b / (1000 * 1024 * 1024) << " GB]" << std::endl;

    if (mpi_rank_ == 0) std::cout << "done" << std::endl;

    // std::cout << ecc_sh_ptr_->print_atmos_parameters(0,0,1);
    // std::cout << ecc_sh_ptr_->print_atmos_parameters(1,1,65);
}


// emissivity from scattering (line + continuum)
void MF_context::update_emission(const Vec &I_vec, const bool approx){ 	

	PetscErrorCode ierr; 
	       
    const auto g_dev   = RT_problem_->space_grid_->view_device();  
    const auto eta_dev = RT_problem_->eta_field_ ->view_device();
    const auto S_dev   = RT_problem_->S_field_   ->view_device(); 

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
    ierr = VecGetOwnershipRange(I_vec, &istart, &iend);CHKERRV(ierr);	
	
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

        // get I field 
        ierr = VecGetValues(I_vec, block_size, ix, &input[0]);CHKERRV(ierr);   

        // compute grid indeces from Vec index i_vec
        i = i_start + counter_i;
        j = j_start + counter_j;
        k = k_start + counter_k;     

        if (i >= i_end) std::cout << "ERROR with counters in update_emission(), i = " << i << std::endl;
        if (j >= j_end) std::cout << "ERROR with counters in update_emission(), j = " << j << std::endl;
        if (k >= k_end) std::cout << "ERROR with counters in update_emission(), k = " << k << std::endl;

    	// set input field
        ecc_sh_ptr_->update_incoming_field(i, j, k, offset_fun_, input.data());
       	
    	if (approx)
    	{
    		const auto out_field = epsilon_fun_approx_(i,j,k);       
	    	rii_include::make_indices_convertion_function<double>(out_field, offset_fun_)(output.data());    
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

        if (counter_i == i_end - g_dev.margin[0])
        {
            counter_i = 0;
            counter_j++;
        }

        if (counter_j == j_end - g_dev.margin[1])
        {
            counter_j = 0;
            counter_k++;
        }
    }    
}


void MF_context::init_serial_fields(const size_t n_tiles){
    
    auto block_size = RT_problem_->block_size_;

    auto N_x = RT_problem_->N_x_;
    auto N_y = RT_problem_->N_y_;
    auto N_z = RT_problem_->N_z_;

    // set the number of local rays and tiles
    n_tiles_ = n_tiles;    
    n_local_rays_ = block_size/mpi_size_;

    if (n_local_rays_ < 4) // TODO: to test
    {
        if (mpi_rank_ == 0) std::cout << "WARNING: mpi_size > number of rays, not tested" << std::endl;
        n_local_rays_ = 4; 
    } 
    else
    {
        if (n_local_rays_ * mpi_size_ != block_size) std::cout << "ERROR in init_serial_fields(): block_size/mpi_size_ not integer" << std::endl;  
    }
    
    
    if (n_local_rays_ % 4 != 0) std::cout << "ERROR in init_serial_fields(): n_local_rays_ should be divisible by 4" << std::endl;        

    tile_size_ = n_local_rays_/n_tiles_;

    if (tile_size_ * n_tiles_ != n_local_rays_) std::cout << "ERROR in init_serial_fields(): n_local_rays_/n_tiles_ not integer" << std::endl;        
    if (tile_size_ % 4 != 0)                    std::cout << "ERROR in init_serial_fields(): tile_size_ should be divisible by 4" << std::endl;            

    // init serial grid
    const bool use_ghost_layers = false;
    if (mpi_rank_ == 0 and use_ghost_layers) std::cout << "\nusing ghost layers for serial grid" << std::endl;    

    space_grid_serial_ = std::make_shared<Grid_t>();    
    space_grid_serial_->init(MPI_COMM_SELF, {(int)N_x, (int)N_y, (int)N_z}, {1, 1, 0}, {}, use_ghost_layers); 

    // create serial fields 
    I_field_serial_   = std::make_shared<Field_t>("I_serial", space_grid_serial_, tile_size_); 
    S_field_serial_   = std::make_shared<Field_t>("S_serial", space_grid_serial_, tile_size_);

    if (n_tiles_ != 1) std::cout << "ERROR: n_tiles_ should be 1 for now (for b indexing eta and rho) " << std::endl;

    eta_field_serial_ = std::make_shared<Field_t>("eta_serial", space_grid_serial_, n_local_rays_); // here could tiles also be used to reduce mem footprint
    rho_field_serial_ = std::make_shared<Field_t>("rho_serial", space_grid_serial_, n_local_rays_);

    // allocate
    I_field_serial_  ->allocate_on_device();     
    S_field_serial_  ->allocate_on_device();     

    eta_field_serial_->allocate_on_device();     
    rho_field_serial_->allocate_on_device();     

    // init remaps 
    I_remap_.init(*(RT_problem_->I_field_), *I_field_serial_);
    S_remap_.init(*(RT_problem_->S_field_), *S_field_serial_);
            
    sgrid::ReMap<Field_t> tmp_remap;

    tmp_remap.init(*(RT_problem_->eta_field_), *eta_field_serial_);
    tmp_remap.from_pgrid_to_pblock(*(RT_problem_->eta_field_), *eta_field_serial_, 0); 

    tmp_remap.init(*(RT_problem_->rho_field_), *rho_field_serial_);
    tmp_remap.from_pgrid_to_pblock(*(RT_problem_->rho_field_), *rho_field_serial_, 0); 
} 


void RT_solver::print_info(){

    // print some output
    if (mpi_rank_ == 0)
    {
        if (mf_ctx_.use_log_interpolation_)
        {
            std::cout << "Using logarithmic interpolation and ";
        }
        else
        {
            std::cout << "Using linear interpolation and ";
        }

        if (mf_ctx_.use_single_long_step_)
        {
            std::cout << "a single step for long rays." << std::endl;
        }
        else
        {
            std::cout << "multiple steps for long rays." << std::endl;
        }

        if (mf_ctx_.use_always_long_ray_) std::cout << "Only long rays are used." << std::endl;

        std::cout << "\n========= Serial formal solver parameters =========" << std::endl;
        std::cout << "n_local_rays = " << mf_ctx_.n_local_rays_ << " (block_size/mpi_size)" << std::endl;             
        std::cout << "tile_size    = " << mf_ctx_.tile_size_    << " (n_local_rays/n_tiles)" << std::endl;    
        std::cout << "n_tiles      = " << mf_ctx_.n_tiles_ << std::endl;                          
        std::cout << "=====================================================" << std::endl;
    } 
}


void RT_solver::assemble_rhs(){ 

    // with test = true data structures are created but not filled
    const bool test = false;

  	if (mpi_rank_ == 0) std::cout           << "\n++++++ Assembling right hand side...+++++++++";
    if (mpi_rank_ == 0 and test ) std::cout << "\n+++++++++++ RHS TEST RHS TEST +++++++++++++";
 
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

    if (not test)
    {     
    	// create rhs field (temporary)
    	auto rhs_field    = std::make_shared<Field_t>("RHS", space_grid, block_size);
    	rhs_field->allocate_on_device(); 	
    	
    	// create eps_th field (temporary)
    	auto eps_th_field = std::make_shared<Field_t>("EPS_TH", space_grid, block_size);
    	eps_th_field->allocate_on_device(); 
    	const auto eps_th_dev = eps_th_field->view_device();	 
    	
    	// fill eps_th =  eps_c_th +  eps_l_th
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

    	// fill rhs_ from formal solve with bc (I_field_serial_ and S_field_serial_ are used as temporary containers here)           	
        mf_ctx_.formal_solve_global(rhs_field, mf_ctx_.I_field_serial_, eps_th_field, mf_ctx_.S_field_serial_, 1.0);       
    	mf_ctx_.field_to_vec(rhs_field, rhs_);  

        // rhs_field->write("rhs_field.raw");   
    }

	if (mpi_rank_ == 0) std::cout << "+++++++++++++++++++++++++++++++++++++++++++++\n"; 	
}


// matrix-free matrix vector multiplication y = (Id - LJ)x
PetscErrorCode UserMult(Mat mat, Vec x, Vec y){

	PetscErrorCode ierr; 

	void *ptr;
   	MatShellGetContext(mat,&ptr);
  	MF_context *mf_ctx_ = (MF_context *)ptr;

  	auto RT_problem = mf_ctx_->RT_problem_;

    Real start = MPI_Wtime();       

    // compute new emission in S_field_ 
    mf_ctx_->update_emission(x);  

    if (RT_problem->mpi_rank_ == 0) printf("update_emission:\t\t%g seconds\n", MPI_Wtime() - start);              
    start = MPI_Wtime();      
  	    
  	// fill rhs_ from formal solve with zero bc  	
	mf_ctx_->formal_solve_global(RT_problem->I_field_, mf_ctx_->I_field_serial_, RT_problem->S_field_, mf_ctx_->S_field_serial_, 0.0);
      
    if (RT_problem->mpi_rank_ == 0) printf("formal_solve_global:\t\t%g seconds\n", MPI_Wtime() - start);              
    
	// copy intensity into petscvec format
	mf_ctx_->field_to_vec(RT_problem->I_field_, y);

	// update I_out = I_in - I_fs (y = x - y)
	ierr = VecAYPX(y, -1.0, x);CHKERRQ(ierr);

  	return ierr;
}


// matrix-free matrix vector multiplication y = (Id - LJ)x
PetscErrorCode UserMult_approx(Mat mat, Vec x, Vec y){

	PetscErrorCode ierr; 

	void *ptr;
   	MatShellGetContext(mat,&ptr);
  	MF_context *mf_ctx_ = (MF_context *)ptr; 

    auto RT_problem = mf_ctx_->RT_problem_; 

    Real start = MPI_Wtime();       

    // compute new emission in S_field_ 
    mf_ctx_->update_emission(x, true);  
   
    if (RT_problem->mpi_rank_ == 0) printf("update CRD emission:\t\t%g seconds\n", MPI_Wtime() - start);              

    // fill rhs_ from formal solve with zero bc     
    mf_ctx_->formal_solve_global(RT_problem->I_field_, mf_ctx_->I_field_serial_, RT_problem->S_field_, mf_ctx_->S_field_serial_, 0.0);
    
    // copy intensity into petscvec format
    mf_ctx_->field_to_vec(RT_problem->I_field_, y);
    
    // update I_out = I_in - I_fs (y = x - y)
    ierr = VecAYPX(y, -1.0, x);CHKERRQ(ierr);

  	return ierr;
}


PetscErrorCode MF_pc_Destroy(PC pc){

	PetscErrorCode ierr;

	MF_context *mf_ctx;

	ierr = PCShellGetContext(pc,(void**)&mf_ctx); CHKERRQ(ierr);   
    ierr = PetscFree(mf_ctx);

    // TODO destroy?

	return ierr;
}

PetscErrorCode MF_pc_Apply(PC pc,Vec x,Vec y){

	PetscErrorCode ierr;

	MF_context *mf_ctx;

	ierr = PCShellGetContext(pc,(void**)&mf_ctx);CHKERRQ(ierr);   

	// apply	
	ierr = KSPSolve(mf_ctx->pc_solver_, x, y);CHKERRQ(ierr);

    // print  iterations 
    int iterations;
    ierr = KSPGetIterationNumber(mf_ctx->pc_solver_, &iterations);CHKERRQ(ierr);
    if (mf_ctx->mpi_rank_ == 0) std::cout << "Preconditioner iterations: " << iterations << std::endl;

	return ierr;
}
