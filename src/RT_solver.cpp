#include "RT_solver.hpp"
#include "cpu_clock.h"
#include <string>

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

	PetscInt istart, row;

	double value;
	
	ierr = VecGetOwnershipRange(v, &istart, NULL);CHKERRV(ierr);	

	int counter = 0;

	for (int k = k_start; k < k_end; ++k)					
	{															
		for (int j = j_start; j < j_end; ++j)
		{
			for (int i = i_start; i < i_end; ++i)				
			{
				for (int b = 0; b < block_size; b++) 
				{
					// set row index and corresponding entry
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
	if (mpi_rank_ == 0) std::cout << "\nCopying Vec to field...";

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

	PetscInt istart, row;

	double value;
	
	ierr = VecGetOwnershipRange(v, &istart, NULL);CHKERRV(ierr);	

	int counter = 0;

	for (int k = k_start; k < k_end; ++k)					
	{															
		for (int j = j_start; j < j_end; ++j)
		{
			for (int i = i_start; i < i_end; ++i)				
			{
				for (int b = 0; b < block_size; b++) 
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

    if (mpi_rank_ == 0) std::cout << "done" << std::endl;
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
            
            for (int b = 0; b < block_size; b = b + 4) 
            {
                I_field_dev.block(i,j,k)[b] = W_T_deep;                
            }            
        }
    }); 
}


void MF_context::find_intersection(double theta, double chi, const double Z_down, const double Z_top, const double L, t_intersect *T) 
{
    // check theta and possibly correct
    if (theta == 0 or chi == 0) 
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
    	if (T->w[i] < -1e-15 or T->w[i] > 1)
    	{
    		std::cout << "WARNING in find_intersection(): w has a problem!" << std::endl;         	
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
            std::cout << "WARNING in find_prolongation(): w has a problem!" << std::endl;          
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


// given a intersection type with N cells and grid indeces ijk, get I1, S1, K1 i.e. quantities needed for the last step of formal solution
std::vector<double> MF_context::long_ray_steps(const std::vector<t_intersect> T, 
                                               const Field_ptr_t I_field, const Field_ptr_t S_field, 
                                               const int i, const int j, const int k, const int block_index)
{   
    const auto N_x = RT_problem_->N_x_;
    const auto N_y = RT_problem_->N_y_;
    
    const auto eta_dev = (formal_solution_Omega_) ? eta_field_serial_Omega_->view_device() : eta_field_serial_->view_device(); 
    const auto rho_dev = (formal_solution_Omega_) ? rho_field_serial_Omega_->view_device() : rho_field_serial_->view_device(); 
    
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

            weight = 0;
            for (int face_vertices = 0; face_vertices < 4; ++face_vertices)
            {
                weight += T[cell].w[face_vertices];
            }            

			for (int face_vertices = 0; face_vertices < 4; ++face_vertices)
			{
				i_intersect = i + T[cell].ix[face_vertices];
				j_intersect = j + T[cell].iy[face_vertices];
				k_intersect = k - T[cell].iz[face_vertices]; 

			    // correction for periodic BC 
                i_intersect = apply_periodic_bc(i_intersect, N_x);
                j_intersect = apply_periodic_bc(j_intersect, N_y);                

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
                i_intersect = apply_periodic_bc(i_intersect, N_x);
                j_intersect = apply_periodic_bc(j_intersect, N_y);                

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

        // compute current interval distance
        const double cell_distance = (cell < N - 1) ? T[cell].distance - T[cell + 1].distance : T[cell].distance;         

		// optical depth step		        
		dtau = coeff * (eta_I_1 + etas[0]) * cell_distance;									

		if (dtau > 0)  std::cout << "ERROR in dtau sign, dtau = " << dtau << std::endl;
        if (dtau == 0) std::cout << "WARNING: dtau = 0, possible e.g. for N_chi = 4"<< std::endl;
        
		formal_solver_.one_step(dtau, K1, K2, S1, S2, I1, I2);

        // // test
        // // get indeces
        // std::vector<int> local_idx;
        // local_idx = RT_problem_->block_to_local(tile_size_* mpi_rank_ + block_index);
        
        // const int j_theta = local_idx[0];
        // const int k_chi   = local_idx[1];
        // const int n_nu    = local_idx[2];

        // const auto mu_grid    = RT_problem_->mu_grid_;
        // const auto theta_grid = RT_problem_->theta_grid_;   
        // const auto chi_grid   = RT_problem_->chi_grid_;   
    
        // const double theta = theta_grid[j_theta];
        // const double mu    = mu_grid[j_theta];     
        // const double chi   = chi_grid[k_chi];                       

        // std::cout << "k = " << k << std::endl;    
        // std::cout << "j_theta = " << j_theta << std::endl;
        // std::cout << "k_chi = "   << k_chi << std::endl;
        // std::cout << "n_nu = "    << n_nu << std::endl;

        //if (j_theta == 4 and k_chi == 0 and n_nu == 45)
        //{                                                                            
            // if (long_ray) std::cout << "WARNING LONG RAY: look at long ray routines for data!" << std::endl;                                               
            // std::cout << "\nk = " << k << std::endl;                                                          
            // std::cout << "cell = "  << cell << std::endl;
            // std::cout << "mu = " << mu << std::endl;
            // std::cout << "chi = "   << chi << std::endl;
            // std::cout << "n_nu = " << n_nu << std::endl;
            // std::cout << "mu = " << mu << std::endl;
            // std::cout << "n = "  << n << std::endl;                                            
            // std::cout << "dz = "<< dz << std::endl;
            
            // std::cout << "mpi_rank_ = " << mpi_rank_ << std::endl;   
            // std::cout << "k_global = " << g_dev.global_coord(2, k_aux) << std::endl;                                                          

            // std::cout << "I1 = "   << I1[0] << std::endl;   
            // std::cout << "Q1 = "   << I1[1] << std::endl;   
            // std::cout << "U1 = "   << I1[2] << std::endl;   
            // std::cout << "V1 = "   << I1[3] << std::endl;   

            // // std::cout << "I2 = "   << I2[0] << std::endl;    
            // // std::cout << "Q2 = "   << I2[1] << std::endl;   
            // std::cout << "U2 = "   << I2[2] << std::endl;   
            // std::cout << "V2 = "   << I2[3] << std::endl;   

            // std::cout << "S1[2] = "   << S1[2] << std::endl;   
            // std::cout << "S2[2] = "   << S2[2] << std::endl;   

            // std::cout << "S1[3] = "   << S1[3] << std::endl;   
            // std::cout << "S2[3] = "   << S2[3] << std::endl;   

            // std::cout << "S1 = " << std::endl;
            // for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S1[i_stokes] << std::endl;

            // std::cout << "S2 = " << std::endl;
            // for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S2[i_stokes] << std::endl;

            // std::cout << "K1 = " << std::endl;
            // for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K1[i_stokes] << std::endl;

            // std::cout << "K2 = " << std::endl;
            // for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K2[i_stokes] << std::endl;

            // std::cout << "dtau = " << dtau  << std::endl;    
            // std::cout << "etas[0] = " << etas[0] << std::endl;  
            // std::cout << "eta_I_1 = " << eta_I_1 << std::endl;                     
            // std::cout << "distance = " << cell_distance << std::endl; 
        //}       
	}	
                                                                                                                    
    return I2;
}


// given a intersection type with N cells and grid indeces ijk, get I1, S1, K1 i.e. quantities needed for the last step of formal solution
std::vector<double> MF_context::long_ray_steps_quadratic(const std::vector<t_intersect> T, 
                                                         const Field_ptr_t I_field, const Field_ptr_t S_field, 
                                                         const int i, const int j, const int k, const int block_index,
                                                         bool print_flag) // to test
{                         
    if (not use_always_long_ray_) std::cout << "WARNING: short ray in long_ray_steps_quadratic()!" << std::endl;    

    const auto N_x = RT_problem_->N_x_;
    const auto N_y = RT_problem_->N_y_;

    const auto eta_dev = (formal_solution_Omega_) ? eta_field_serial_Omega_->view_device() : eta_field_serial_->view_device(); 
    const auto rho_dev = (formal_solution_Omega_) ? rho_field_serial_Omega_->view_device() : rho_field_serial_->view_device(); 
    
    const auto I_dev = I_field->view_device();     
    const auto S_dev = S_field->view_device(); 
    
    // coeff trap + cm conversion = - 0.5 * 1e5;
    const double coeff = -50000;

    // number of traversec cells 
    const int N = T.size() - 1;

    int i_intersect, j_intersect, k_intersect, b_index;

    double eta_I_1, weight, dtau_1, dtau_2, cell_distance;

    double distance_test, etas_1_print;

    std::vector<double> I1(4), I2(4), S1(4), S2(4), S3(4), etas(4), rhos(4), K1(16), K2(16), K3(16);     
   
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
                i_intersect = apply_periodic_bc(i_intersect, N_x);
                j_intersect = apply_periodic_bc(j_intersect, N_y);
                
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

            // save for later use
            eta_I_1 = etas[0];

            if (print_flag) etas_1_print =  etas[0];
        }
        else // reuse quantities in (2) 
        {
            S1 = S2;
            I1 = I2;
            K1 = K2;        
        }        

        ////////////////////////////////////////////////////////////////////////////

        // quantities in (2) 
        if (cell == 0) 
        {            
            if (cell == N - 1) // in case of short ray this is also the last iterate
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

                // compute current interval distance (different formula for last cell)
                cell_distance = T[cell].distance;
            }
            else 
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
                    i_intersect = apply_periodic_bc(i_intersect, N_x);
                    j_intersect = apply_periodic_bc(j_intersect, N_y);                

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

                // compute current interval distance
                cell_distance = T[cell].distance - T[cell + 1].distance;
            }

            K2 = assemble_propagation_matrix_scaled(etas, rhos);     
            
            dtau_1 = coeff * (eta_I_1 + etas[0]) * cell_distance;                  

            // distance_test = cell_distance;
            
            if (dtau_1 > 0)  std::cout << "ERROR in dtau_1 sign, dtau_1 = " << dtau_1 << std::endl;   
            if (dtau_1 == 0) std::cout << "WARNING: dtau_1 = 0, possible e.g. for N_chi = 4" << std::endl;                                                         
        }
        else // reuse
        {
            S2 = S3;            
            K2 = K3;

            dtau_1 = dtau_2;
        }

        // for next tau
        eta_I_1 = etas[0];

        ////////////////////////////////////////////////////////////////////////////

        // quantities in (3)
        if (cell == N - 2) // no interpolation
        {
            for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
            {
                b_index = block_index + i_stokes;                                       
        
                // get eta and rho
                etas[i_stokes] = eta_dev.block(i,j,k)[b_index]; 
                rhos[i_stokes] = rho_dev.block(i,j,k)[b_index];

                // set S3
                S3[i_stokes] = S_dev.block(i,j,k)[b_index];                                                         
            }     

            // compute current interval distance 
            cell_distance = T[cell + 1].distance;          
        }
        else
        {
            // init
            for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
            {
                etas[i_stokes] = 0;
                rhos[i_stokes] = 0;
                S3[i_stokes]   = 0;            
            }

            for (int face_vertices = 0; face_vertices < 4; ++face_vertices)
            {
                const int next_cell = (cell == N - 1) ? cell + 1 : cell + 2; 

                i_intersect = i + T[next_cell].ix[face_vertices];
                j_intersect = j + T[next_cell].iy[face_vertices];
                k_intersect = k - T[next_cell].iz[face_vertices]; 

                // correction for periodic boundary
                i_intersect = apply_periodic_bc(i_intersect, N_x);
                j_intersect = apply_periodic_bc(j_intersect, N_y);                

                weight = T[next_cell].w[face_vertices];  

                for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
                {
                    b_index = block_index + i_stokes;                                       

                    // get eta and rho
                    etas[i_stokes] += weight * eta_dev.block(i_intersect,j_intersect,k_intersect)[b_index]; 
                    rhos[i_stokes] += weight * rho_dev.block(i_intersect,j_intersect,k_intersect)[b_index];                    

                    // set S3
                    S3[i_stokes] += weight * S_dev.block(i_intersect,j_intersect,k_intersect)[b_index];                                                         
                }               
            }   

            // compute current interval distance
            if (cell < N - 1)
            {
                cell_distance = T[cell + 1].distance - T[cell + 2].distance;         
            } 
            else //different mechanism in last cell
            {
                cell_distance = T[cell + 1].distance;         
            }            
        }

        K3 = assemble_propagation_matrix_scaled(etas, rhos);
                                
        // optical depth step               
        dtau_2 = coeff * (eta_I_1 + etas[0]) * cell_distance; 
       
        if (dtau_2 > 0)  std::cout << "ERROR in dtau_2 sign, dtau_2 = " << dtau_2 << std::endl;  
        if (dtau_2 == 0) std::cout << "WARNING: dtau2 = 0, possible e.g. for N_chi = 4" << std::endl;
        
        formal_solver_.one_step_quadratic(dtau_1, dtau_2, K1, K2, K3, S1, S2, S3, I1, I2);       

        // TEST
        bool print_flag2 = true;

        // Real mu, chi;

        // if (not formal_solution_Omega_)
        // {            
        //     // get indeces
        //     std::vector<int> local_idx;
        //     local_idx = RT_problem_->block_to_local(tile_size_* mpi_rank_ + block_index);
            
        //     const int j_theta = local_idx[0];
        //     const int k_chi   = local_idx[1];            
        //     const int n_nu    = local_idx[2];

        //     if (j_theta == 7 and k_chi == 15 and n_nu == 20)
        //     {
        //         print_flag2 = true;

        //         const auto mu_grid  = RT_problem_->mu_grid_;        
        //         const auto chi_grid = RT_problem_->chi_grid_; 
                    
        //         mu  = mu_grid[j_theta];     
        //         chi = chi_grid[k_chi];  
        //     }            
        // }
        // else
        // {
        //     print_flag2 = true;
        // }        
        
        // const auto N_theta = RT_problem_->N_theta_;               

        // bool nu_n_equal_zero = (block_index == tile_size_/2);
        
        if (print_flag and print_flag2)            
        {
            std::cout << "mpi_rank_ = " << mpi_rank_ << std::endl;   
        //  std::cout << "j_theta = " << j_theta << std::endl;
        //  std::cout << "k_chi = "   << k_chi << std::endl;         

            // // }
        

        // if (j_theta == N_theta - 1 and k_chi == 0 and n_nu == 0 and i == 0 and j == 0)
        // {                                                                            
            // if (not formal_solution_Omega_)
            // {
            //     // std::cout << "theta = " << theta << std::endl;
            //     std::cout << "mu = "   << mu << std::endl;
            //     std::cout << "chi = "  << chi << std::endl;            
            // }
        //     // std::cout << "n = "  << n << std::endl;                                            
        //     // std::cout << "dz = "<< dz << std::endl;           
           
            // std::cout << "mpi_rank_ = " << mpi_rank_ << std::endl;   
            std::cout << "k = " << k << std::endl;                                              
            std::cout << "block_index = " << block_index << std::endl;                                              
            

                // const auto T_dev = RT_problem_->T_->view_device();

                // std::cout << "i_in = "   << i  << std::endl;       
                // std::cout << "j_in = "   << j  << std::endl;       
                // std::cout << "k_in = "   << k  << std::endl;       

                // std::cout << "T = "   << T_dev.ref(i,j,k)  << std::endl;                   
                
                std::cout << "I1 = "   << I1[0] << std::endl;   
                // std::cout << "Q1 = "   << I1[1] << std::endl;   
                // std::cout << "U1 = "   << I1[2] << std::endl;   
                // std::cout << "V1 = "   << I1[3] << std::endl;   

                std::cout << "I2 = "   << I2[0] << std::endl;    
                // std::cout << "Q2 = "   << I2[1] << std::endl;   
                // std::cout << "U2 = "   << I2[2] << std::endl;   
                // std::cout << "V2 = "   << I2[3] << std::endl; 

                // std::cout << "dtau_1 = " << dtau_1  << std::endl;    
                // std::cout << "dtau_2 = " << dtau_2  << std::endl;   

                // std::cout << "distance_1 = " << distance_test << std::endl;
                // std::cout << "distance_2 = " << cell_distance << std::endl;             

                std::cout << "S1[0] = " <<  S1[0] << std::endl;
                std::cout << "S2[0] = " <<  S2[0] << std::endl;
                std::cout << "S3[0] = " <<  S3[0] << std::endl;            
                std::cout << "K1[0] = " <<  K1[0] << std::endl;
                std::cout << "K2[0] = " <<  K2[0] << std::endl;
                std::cout << "K3[0] = " <<  K3[0] << std::endl;            


            // std::cout << "S1 = " << std::endl;
            // for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S1[i_stokes] << std::endl;

            // std::cout << "S2 = " << std::endl;
            // for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S2[i_stokes] << std::endl;

            // std::cout << "S3 = " << std::endl;
            // for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S3[i_stokes] << std::endl;
                                            
            // std::cout << "K1 = " << std::endl;
            // for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K1[i_stokes] << std::endl;

            // std::cout << "K2 = " << std::endl;
            // for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K2[i_stokes] << std::endl;   

            // std::cout << "K3 = " << std::endl;
            // for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K3[i_stokes] << std::endl;   
            
            // std::cout << "etas_1 = " << etas_1_print << std::endl;  
            // std::cout << "etas_2 = " << eta_I_1 << std::endl;  
            // std::cout << "etas_3 = " << etas[0] << std::endl;                                 
        }       
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

    const auto eta_dev = (formal_solution_Omega_) ? eta_field_serial_Omega_->view_device() : eta_field_serial_->view_device(); 
    const auto rho_dev = (formal_solution_Omega_) ? rho_field_serial_Omega_->view_device() : rho_field_serial_->view_device(); 

    const auto I_dev = I_field->view_device();     
    const auto S_dev = S_field->view_device(); 

    // coeff trap + cm conversion = - 0.5 * 1e5;
    const double coeff = -50000;
    
    int i_intersect, j_intersect, k_intersect, b_index;

    double eta_I_1, weight;
    double total_distance = 0;

    std::vector<double> I1(4), I2(4), S1(4), S2(4), etas(4), rhos(4), K1(16), K2(16);
   
    // total distance 
    total_distance = T[0].distance;  

    // quantities in (1)  
    for (int i_stokes = 0; i_stokes < 4; ++i_stokes) 
    {
        // interpolate
        etas[i_stokes] = 0;
        rhos[i_stokes] = 0;
        S1[i_stokes]   = 0;
        I1[i_stokes]   = 0;
    }

    const double debug_value = std::abs(T[0].iz[0] + T[0].iz[1] + T[0].iz[2] + T[0].iz[3]);

    if (debug_value != 4) std::cout << "ERROR in single_long_ray_step()" << std::endl;

    for (int face_vertices = 0; face_vertices < 4; ++face_vertices)
    {

        i_intersect = i + T[0].ix[face_vertices];
        j_intersect = j + T[0].iy[face_vertices];
        k_intersect = k - T[0].iz[face_vertices]; 
        
        // correction for periodic BC 
        i_intersect = apply_periodic_bc(i_intersect, N_x);
        j_intersect = apply_periodic_bc(j_intersect, N_y);                
       
        weight = T[0].w[face_vertices];          
        
        for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
        {
            b_index = block_index + i_stokes;                                       
        
            // interpolate eta and rho
            etas[i_stokes] += weight * eta_dev.block(i_intersect,j_intersect,k_intersect)[b_index]; 
            rhos[i_stokes] += weight * rho_dev.block(i_intersect,j_intersect,k_intersect)[b_index];

            // interpolate S1 and I1
            S1[i_stokes] += weight * S_dev.block(i_intersect,j_intersect,k_intersect)[b_index];                                             
            I1[i_stokes] += weight * I_dev.block(i_intersect,j_intersect,k_intersect)[b_index];                
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

    if (dtau > 0 ) std::cout << "ERROR in dtau sign, dtau = " << dtau << std::endl;
    if (dtau == 0) std::cout << "WARNING: dtau = 0, possible e.g. for N_chi = 4"<< std::endl;
    
    formal_solver_.one_step(dtau, K1, K2, S1, S2, I1, I2);    
                                                                                                                            
    return I2;
}


void MF_context::formal_solve_ray(const Real theta, const Real chi)
{       
    const Real mu = cos(theta);

    if (mpi_rank_ == 0) std::cout << "\nStart formal solution for mu = " << mu << 
                                    ", theta = " << theta << ", and chi = " << chi << std::endl;    

    // init some quantities         
    const auto N_x = RT_problem_->N_x_;
    const auto N_y = RT_problem_->N_y_;
    const auto N_z = RT_problem_->N_z_;
    
    const auto N_nu = RT_problem_->N_nu_;

    const auto block_size = 4 * N_nu;
    
    const auto depth_grid = RT_problem_->depth_grid_;   
    const auto L          = RT_problem_->L_;            

    const auto eta_dev = eta_field_serial_Omega_->view_device(); 
    const auto rho_dev = rho_field_serial_Omega_->view_device(); 

    const auto I_dev = I_field_serial_Omega_->view_device();      
    const auto S_dev = S_field_serial_Omega_->view_device();      

    const auto g_dev = space_grid_serial_->view_device();   

    // indeces
    const int i_start = g_dev.margin[0]; 
    const int j_start = g_dev.margin[1];
    const int k_start = g_dev.margin[2];

    if (i_start > 0 or j_start > 0 or k_start > 0) std::cout << "WARNING: periodic BC hardcoded for margin = 0!" << std::endl;

    const int i_end = i_start + g_dev.dim[0];
    const int j_end = j_start + g_dev.dim[1];
    const int k_end = k_start + g_dev.dim[2];           

    const int stencil_size = formal_solver_.get_stencil_size();

    bool use_linear_method;
    
    int i_aux, j_aux, k_aux, k_global, b_index;

    std::vector<int> i_intersect(4), j_intersect(4), k_intersect(4);

    // serial indexing coeffs    
    std::vector<int> local_idx;    

    // misc coeffs
    double dtau, weight, eta_I_1, dz;
    
    bool boundary, horizontal_face, long_ray;

    // quantities depending on spatial point i
    std::vector<double> I1(4), I2(4), S1(4), S2(4), etas(4), rhos(4), K1(16), K2(16);

    // intersection object
    t_intersect intersection_data;        
    std::vector<t_intersect> intersection_data_long_ray;

    // for next interesection along the outgoing ray, used just in case of stencil_size = 3
    t_intersect intersection_data_next;

    // minus for optical depth conversion, trap rule and conversion to cm (- 0.5 * 1e5)
    const double coeff = -50000;
    
    double comm_timer     = 0;
    double one_step_timer = 0;
    double total_timer    = 0;

    const bool timing_debug = false;

    if (timing_debug) MPI_Barrier(MPI_COMM_WORLD);
    Real start_total = MPI_Wtime(); 

    ///////////// data movement ////////////////////
    
    // write eta and rhos to the serial grid
    sgrid::ReMap<Field_t> tmp_remap;

    tmp_remap.init(*(RT_problem_->eta_field_Omega_), *eta_field_serial_Omega_);
    tmp_remap.from_pgrid_to_pblock(*(RT_problem_->eta_field_Omega_), *eta_field_serial_Omega_, 0); 

    tmp_remap.init(*(RT_problem_->rho_field_Omega_), *rho_field_serial_Omega_);
    tmp_remap.from_pgrid_to_pblock(*(RT_problem_->rho_field_Omega_), *rho_field_serial_Omega_, 0); 
    
    // write S to the serial grid
    S_remap_Omega_.from_pgrid_to_pblock(*RT_problem_->S_field_Omega_, *S_field_serial_Omega_, 0);                                    

    ////////////////////////////////////////////////

    // const bool idle_proc = (mpi_rank_ * local_block_size_ > block_size - 1);

    // TODO FIX
    const bool idle_proc = (eta_dev.block(i_start,j_start,k_start)[0] == 0);
                         
    // communication timer                 
    if (timing_debug) MPI_Barrier(MPI_COMM_WORLD);
    Real start_comm = MPI_Wtime();                                    
            
    comm_timer += MPI_Wtime() - start_comm;          
    
    if (not idle_proc)
    {                       
        // /////////// TEST //////////////            
        // if (g_dev.global_coord(0, i_start) == 0 and
        //     g_dev.global_coord(1, j_start) == 0 and 
        //     g_dev.global_coord(2, k_start) == 0)
        // {                    
        //     std::cout << "--------------- mpi_rank_ = " << mpi_rank_ << std::endl;

        //     // for (int bb = 0; bb < local_block_size_; +bb)
        //     // {
        //         std::cout << "eta_serial = " << eta_dev.block(i_start,j_start,k_start)[0] << std::endl;
        //     // }

        //     std::cout << "---------------"<< std::endl;
        // }

        // const auto eta_dev_par = RT_problem_->eta_field_Omega_->view_device();     
        // const auto g_dev_par   = RT_problem_->space_grid_->view_device();   

        // const int i_s = g_dev_par.margin[0];
        // const int j_s = g_dev_par.margin[1];
        // const int k_s = g_dev_par.margin[2];
        
        // if (g_dev_par.global_coord(0, i_s) == 0 and g_dev_par.global_coord(1, j_s) == 0 and g_dev_par.global_coord(2, k_s) == 0)
        // {        
        //     std::cout << "---------------"<< std::endl;

        //     for (int bb = 0; bb < 4*N_nu; bb = bb + 4)
        //     {
        //         std::cout << "eta = " << eta_dev_par.block(i_s,j_s,k_s)[bb] << std::endl;
        //     }

        //     std::cout << "---------------"<< std::endl;
        // }

        ////////////////////


        // loops over spatial points
        for (int k = k_start; k < k_end; ++k)                   
        {                                                  
            for (int j = j_start; j < j_end; ++j)
            {                
                for (int i = i_start; i < i_end; ++i)
                {                          
                    k_aux = (mu > 0.0) ? k_end - 1 - k + g_dev.margin[2]: k; 

                    // depth index
                    k_global = g_dev.global_coord(2, k_aux);                             

                    boundary = (k_global == 0 and mu < 0) or (k_global == N_z - 1 and mu > 0);
                    
                    if (not boundary)
                    {                       
                        // set vertical box size
                        dz = (mu > 0) ? depth_grid[k_global] -  depth_grid[k_global + 1] : depth_grid[k_global - 1] - depth_grid[k_global];                                                                
                            
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

                        // for quadratic stencil consider an extra intersection point (on the boundary linear is used)
                        use_linear_method = (stencil_size == 2);
                        
                        if (stencil_size == 3)                            
                        {                                                 
                            if ( k_global > 0 and k_global < N_z - 1)
                            {
                                const double dz_2    = (mu > 0) ? depth_grid[k_global - 1] -  depth_grid[k_global] : depth_grid[k_global] - depth_grid[k_global + 1]; 
                                const double theta_2 = PI - theta;
                                const double chi_2   = chi + PI;                                           

                                // find one extra intersection data for last stencil point
                                find_intersection(theta_2, chi_2, dz_2, dz_2, L, &intersection_data_next);       

                                // put all the intersection data in one vector       
                                if (not long_ray) 
                                {
                                    intersection_data_long_ray.clear();
                                    intersection_data_long_ray.push_back(intersection_data);
                                }

                                intersection_data_long_ray.push_back(intersection_data_next);                                               
                            }
                            else
                            {
                                // use a linear method in one_step()                                            
                                use_linear_method = true;                                         
                            }
                        }                                    
                        
                        // set intersection indeces 
                        if (not long_ray) // TODO include and use_linear_method
                        {
                            for (int face_v = 0; face_v < 4; ++face_v)
                            {
                                i_intersect[face_v] = i_aux + intersection_data.ix[face_v];
                                j_intersect[face_v] = j_aux + intersection_data.iy[face_v];
                                k_intersect[face_v] = k_aux - intersection_data.iz[face_v]; // minus because k increases going downwards  
                                
                                // impose periodic BC
                                i_intersect[face_v] = apply_periodic_bc(i_intersect[face_v], N_x);
                                j_intersect[face_v] = apply_periodic_bc(j_intersect[face_v], N_y);                                                            
                            }                                                           
                        }                                                                                                             

                        // loop over block (frequencies)
                        for (int b = 0; b < local_block_size_; b = b + 4)
                        {                
                            if (timing_debug) MPI_Barrier(MPI_COMM_WORLD);                                                                                                                 
                            Real start_one = MPI_Wtime();                                               
                            
                            // solve ODE
                            if (not use_linear_method)
                            {
                                if (use_single_long_step_ and mpi_rank_ == 0) std::cerr << "WARNING: use_single_long_step_ not supported with BESSER" << std::endl;

                                // for debug
                                bool print_flag = false;

                                // if (mpi_rank_ == 48 and g_dev.global_coord(0, i_aux) == 0 and g_dev.global_coord(1, j_aux) == 0) print_flag = true;
                           
                                
                                I2 = long_ray_steps_quadratic(intersection_data_long_ray, I_field_serial_Omega_, S_field_serial_Omega_,
                                                                                            i_aux, j_aux, k_aux, b, print_flag);                                            
                            }
                            else if (long_ray)
                            {          
                                if (use_single_long_step_)
                                {
                                    I2 = single_long_ray_step(intersection_data_long_ray, I_field_serial_Omega_, S_field_serial_Omega_,
                                                                                            i_aux, j_aux, k_aux, b);
                                }
                                else
                                {
                                    I2 = long_ray_steps(intersection_data_long_ray, I_field_serial_Omega_, S_field_serial_Omega_,
                                                                                            i_aux, j_aux, k_aux, b);
                                }                                                                         
                            }
                            else // short ray
                            {                                                   
                                // set S2
                                for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
                                {               
                                    b_index = b + i_stokes;

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
                                    // interpolate                              
                                    etas[i_stokes] = 0;
                                    rhos[i_stokes] = 0;
                                    S1[i_stokes]   = 0;
                                    I1[i_stokes]   = 0;                                             
                                }
                            
                                // loop over the four vertex of the intersection face
                                for (int face_v = 0; face_v < 4; ++face_v)
                                {                                                       
                                    weight = intersection_data.w[face_v];
                                
                                    for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
                                    {
                                        b_index = b + i_stokes;                                       

                                        // interpolate eta and rho                                                  
                                        etas[i_stokes] += weight * eta_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index]; 
                                        rhos[i_stokes] += weight * rho_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index];                                                    

                                        // interpolate S1 and I1
                                        S1[i_stokes] += weight * S_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index];                                                        
                                        I1[i_stokes] += weight * I_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index];                                                        
                                    }                                           
                                }                                                                                                                                                   

                                K1 = assemble_propagation_matrix_scaled(etas, rhos);                                                                            
                                
                                // optical depth step                               
                                dtau = coeff * (eta_I_1 + etas[0]) * intersection_data.distance;                                    
                                
                                if (dtau > 0)  std::cout << "ERROR in dtau sign, dtau = " << dtau << std::endl;  
                                if (dtau == 0) std::cout << "WARNING: dtau = 0, possible e.g. for N_chi = 4"<< std::endl;

                                formal_solver_.one_step(dtau, K1, K2, S1, S2, I1, I2);

                                // // test
                                // if (i == i_start and j == j_start and mpi_rank_ == 48)                                                
                                // {                                                                                                                         
                                //     std::cout << "\nk = " << k << std::endl;                                              
                                //     std::cout << "\nk_global = " << g_dev.global_coord(2, k_aux) << std::endl;                                              
                                //     // std::cout << "theta = " << theta << std::endl;
                                //     std::cout << "mu = "  << mu  << std::endl;                                                
                                //     std::cout << "chi = " << chi << std::endl;                                                
                                //     // std::cout << "n = "  << n << std::endl;   
                                //     // std::cout << "b = "  << b << std::endl;                                                                                            
                                    
                                //     // std::cout << "dtau = "     << dtau  << std::endl;   
                                //     // std::cout << "coeff = "    << coeff << std::endl;   
                                //     // std::cout << "etas[0] = "  << etas[0] << std::endl;  
                                //     // std::cout << "eta_I_1 = "  << eta_I_1 << std::endl;                                                   
                                //     // std::cout << "distance = " << intersection_data.distance << std::endl;                                                                                                           

                                //     std::cout << "mpi_rank_ = " << mpi_rank_ << std::endl;                                               

                                //     std::cout << "I1 = "   << I1[0] << std::endl;   
                                //     // std::cout << "Q1 = "   << I1[1] << std::endl;   
                                //     // std::cout << "U1 = "   << I1[2] << std::endl;   
                                //     // std::cout << "V1 = "   << I1[3] << std::endl;   

                                //     std::cout << "I2 = "   << I2[0] << std::endl;    
                                //     // std::cout << "Q2 = "   << I2[1] << std::endl;   
                                //     // std::cout << "U2 = "   << I2[2] << std::endl;                                                   
                                //     // std::cout << "V2 = "   << I2[3] << std::endl;                                                   
                                    
                                //     std::cout << "S1[0] = " << S1[0]<< std::endl;
                                //     std::cout << "S2[0] = " << S2[0]<< std::endl;
                                //     // std::cout << "S1[3] = " << S1[3]<< std::endl;
                                //     // std::cout << "S2[3] = " << S2[3]<< std::endl;
                                //     // // std::cout << "S1 = " << std::endl;
                                //     // for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S1[i_stokes] << std::endl;

                                //     // std::cout << "S2 = " << std::endl;
                                //     // for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S2[i_stokes] << std::endl;
                                                                
                                //     // std::cout << "K1 = " << std::endl;
                                //     // for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K1[i_stokes] << std::endl;

                                //     // std::cout << "K2 = " << std::endl;
                                //     // for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K2[i_stokes] << std::endl;
                                // }
                            
                                one_step_timer += MPI_Wtime() - start_one;                                                                                                                                                      
                            }   
                                                                      
                            // write result
                            for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
                            {                                                           
                                I_dev.block(i_aux,j_aux,k_aux)[b + i_stokes] = I2[i_stokes];                                      
                            }                                 
                        }                                                   
                    }
                }                
            }               
        }      
    }
                  
    if (timing_debug) MPI_Barrier(MPI_COMM_WORLD);
    start_comm = MPI_Wtime();    
    
    I_remap_Omega_.from_pblock_to_pgrid(*I_field_serial_Omega_, *RT_problem_->I_field_Omega_, 0); 

    comm_timer  += MPI_Wtime() - start_comm;     
    total_timer += MPI_Wtime() - start_total;                     
    
    if (mpi_rank_ == 0)
    {
        printf("comm_timer:\t\t%g seconds\n", comm_timer);
        printf("one_step_timer:\t\t%g seconds\n", one_step_timer);                        
        printf("total_timer:\t\t%g seconds\n", total_timer);                        
    }           
}
      

void MF_context::formal_solve_global(Field_ptr_t I_field, const Field_ptr_t S_field, const Real I0)
{
	if (mpi_rank_ == 0) std::cout << "\nStart global formal solution..." << std::endl;
    
	// init some quantities 	    
    const auto N_x = RT_problem_->N_x_;
    const auto N_y = RT_problem_->N_y_;
    const auto N_z = RT_problem_->N_z_;

    const auto N_theta = RT_problem_->N_theta_; // these can be removed (used for testing)
    const auto N_chi   = RT_problem_->N_chi_;
	
	const auto block_size = RT_problem_->block_size_;
	const auto tot_size   = RT_problem_->tot_size_;
	
	const auto mu_grid    = RT_problem_->mu_grid_;
	const auto theta_grid = RT_problem_->theta_grid_;	
	const auto chi_grid   = RT_problem_->chi_grid_;	
	const auto depth_grid = RT_problem_->depth_grid_;	
	const auto L          = RT_problem_->L_;		

    const auto eta_dev = eta_field_serial_->view_device(); 
    const auto rho_dev = rho_field_serial_->view_device(); 

    const auto I_dev = I_field_serial_->view_device();		
	const auto S_dev = S_field_serial_->view_device();	    

	const auto g_dev = space_grid_serial_->view_device();   

	// indeces
	const int i_start = g_dev.margin[0]; 
	const int j_start = g_dev.margin[1];
	const int k_start = g_dev.margin[2];

	if (i_start > 0 or j_start > 0 or k_start > 0) std::cout << "WARNING: periodic BC hardcoded for margin = 0!" << std::endl;

	const int i_end = i_start + g_dev.dim[0];
	const int j_end = j_start + g_dev.dim[1];
	const int k_end = k_start + g_dev.dim[2];	    

    const int stencil_size = formal_solver_.get_stencil_size();

    bool use_linear_method;
    
	int i_aux, j_aux, k_aux, k_global, b_start, b_index;

    std::vector<int> i_intersect(4), j_intersect(4), k_intersect(4);

    // serial indexing coeffs    
    std::vector<int> local_idx;
    int block_start, block_end, j_theta_start, k_chi_start, n_nu_start, j_theta_end, k_chi_end, n_nu_end;

	// misc coeffs
	double theta, chi, mu, dtau, weight, eta_I_1, dz;
	
	bool boundary, horizontal_face, long_ray;

	// quantities depending on spatial point i
	std::vector<double> I1(4), I2(4), S1(4), S2(4), etas(4), rhos(4), K1(16), K2(16);

	// intersection object
	t_intersect intersection_data;        
    std::vector<t_intersect> intersection_data_long_ray;

    // for next interesection along the outgoing ray, used just in case of stencil_size = 3
    t_intersect intersection_data_next;

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
        
        // write S to the serial grid and I to get initial condition              
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

            if (j_theta_end < j_theta_start) { std::cout << "ERROR with j_theta partition: N_theta*(N_chi)*[N_dirs] shoud be divisible by mpi_size (using extra parentesis ()[] as mpi_size increases)!" << std::endl; throw_error = true; }      
            if (k_chi_end   < k_chi_start)   { std::cout << "ERROR with k_chi partition: N_theta*(N_chi)*[N_dirs] shoud be divisible by mpi_size! (using extra parentesis ()[] as mpi_size increases)!"  << std::endl; throw_error = true; }     
            if (n_nu_end    < n_nu_start)    { std::cout << "ERROR with n_nu partition: N_theta*(N_chi)*[N_dirs] shoud be divisible by mpi_size! (using extra parentesis ()[] as mpi_size increases)!"   << std::endl; throw_error = true; }   
        
            if (local_idx[3] != 3) { std::cout << "ERROR in block decomposition in formal_solve_global(), i_stokes_end not 3" << std::endl; throw_error = true; }      
            
            if (throw_error) throw "ERROR with block decomposition";
        
    		// loop over spatial points
    		for (int k = k_start; k < k_end; ++k)					
    		{									            
    			for (int j = j_start; j < j_end; ++j)
    			{
    				for (int i = i_start; i < i_end; ++i)
    				{					                       
    					// loop over directions 
    					for (int j_theta = j_theta_start; j_theta < j_theta_end; ++j_theta)
    					{
    						theta = theta_grid[j_theta];
    						mu    = mu_grid[j_theta];						

    						k_aux = (mu > 0.0) ? k_end - 1 - k + g_dev.margin[2]: k; 

    						// depth index
    						k_global = g_dev.global_coord(2, k_aux);	                         

    						boundary = (k_global == 0 and mu < 0) or (k_global == N_z - 1 and mu > 0);
    						
    						if (not boundary)
    						{						
    							// set vertical box size
    							dz = (mu > 0) ? depth_grid[k_global] -  depth_grid[k_global + 1] : depth_grid[k_global - 1] - depth_grid[k_global];                                                                

    							for (int k_chi = k_chi_start; k_chi < k_chi_end; ++k_chi)
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

                                    // for quadratic stencil consider an extra intersection point (on the boundary linear is used)
                                    use_linear_method = (stencil_size == 2);
                                    
                                    if (stencil_size == 3)
                                    {                                   
                                        // last point has to use linear method     
                                        if (k_global > 0 and k_global < N_z - 1)
                                        {
                                            const double dz_2    = (mu > 0) ? depth_grid[k_global - 1] -  depth_grid[k_global] : depth_grid[k_global] - depth_grid[k_global + 1]; 
                                            const double theta_2 = PI - theta;
                                            const double chi_2   = chi + PI;                                           

                                            // find one extra intersection data for last stencil point
                                            find_intersection(theta_2, chi_2, dz_2, dz_2, L, &intersection_data_next);       

                                            // put all the intersection data in one vector       
                                            if (not long_ray) 
                                            {
                                                intersection_data_long_ray.clear();
                                                intersection_data_long_ray.push_back(intersection_data);
                                            }

                                            intersection_data_long_ray.push_back(intersection_data_next);                                               
                                        }
                                        else
                                        {
                                            // use a linear method in one_step()                                            
                                            use_linear_method = true;                                         
                                        }
                                    }                                    
                                    
                                    // set intersection indeces 
                                    if (not long_ray) // TODO include  and use_linear_method
                                    {
                                        for (int face_v = 0; face_v < 4; ++face_v)
                                        {
                                            i_intersect[face_v] = i_aux + intersection_data.ix[face_v];
                                            j_intersect[face_v] = j_aux + intersection_data.iy[face_v];
                                            k_intersect[face_v] = k_aux - intersection_data.iz[face_v]; // minus because k increases going downwards  
                                            
                                            // impose periodic BC
                                            i_intersect[face_v] = apply_periodic_bc(i_intersect[face_v], N_x);
                                            j_intersect[face_v] = apply_periodic_bc(j_intersect[face_v], N_y);                                                            
                                        }                                                           
                                    }                                                                 
                                    
    								// loop on freqs
                                    for (int n = n_nu_start; n < n_nu_end; ++n)
    								{			     
                                        if (timing_debug) MPI_Barrier(MPI_COMM_WORLD);                                                                                                                 
                                        Real start_one = MPI_Wtime();                                               

    									// block index (corrected for tile size)                                    
    									b_start = RT_problem_->local_to_block(j_theta, k_chi, n) % tile_size_;                                   

    									// solve ODE
                                        if (not use_linear_method)
                                        {
                                            if (use_single_long_step_ and mpi_rank_ == 0) std::cerr << "WARNING: use_single_long_step_ not supported with BESSER" << std::endl;

                                            // for debug
                                            bool print_flag = false;

                                            // if (g_dev.global_coord(0, i_aux) == 0 and g_dev.global_coord(1, j_aux) == 0) print_flag = true;

                                           // if (j_theta ==  N_theta/2 + 1 and k_chi == 0 and n == 0  and i_aux == 0 and j_aux == 0)                                                                                            
                                            //{
                                              //  std::cout << "\nk = " << k << std::endl;
                                               // std::cout << "theta = " << theta << std::endl;
                                                //std::cout << "mu = " << mu << std::endl;          
                                                //std::cout << "chi = "<< chi << std::endl;      
                                                //std::cout << "i = " << g_dev.global_coord(0, i_aux) << std::endl;   
                                                //std::cout << "j = " << g_dev.global_coord(1, j_aux) << std::endl;   
                                                //std::cout << "k = " << g_dev.global_coord(2, k_aux) << std::endl;   

                                                //print_flag = true;                                          
                                            //}
        
                                            // std::cout << "cells = "   << intersection_data_long_ray.size() << std::endl;         
                                            // std::cout << "b_start = " << b_start << std::endl;
                                            // std::cout << "tile_size_ = " << tile_size_ << std::endl;
                                            // // std::cout << "L = "   << L << std::endl;                                          
                                            // std::cout << "dz = "   << dz << std::endl;         
                                            // if (use_linear_method)  std::cout << "LM"  << std::endl;       
                                            // std::cout << "stencil_size = " << stencil_size << std::endl;                                                                  

                                            
                                            I2 = long_ray_steps_quadratic(intersection_data_long_ray, I_field_serial_, S_field_serial_, i_aux, j_aux, k_aux, b_start, print_flag);                                            
                                        }
    									else if (long_ray)
    									{          
                                            if (use_single_long_step_)
                                            {
                                                I2 = single_long_ray_step(intersection_data_long_ray, I_field_serial_, S_field_serial_, i_aux, j_aux, k_aux, b_start);
                                            }
                                            else
                                            {
                                                I2 = long_ray_steps(intersection_data_long_ray, I_field_serial_, S_field_serial_, i_aux, j_aux, k_aux, b_start);
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
                                                // interpolate								
												etas[i_stokes] = 0;
												rhos[i_stokes] = 0;
												S1[i_stokes]   = 0;
												I1[i_stokes]   = 0;    											
    										}
    									
    										// loop over the four vertex of the intersection face
    										for (int face_v = 0; face_v < 4; ++face_v)
    										{			                                            
                                                weight = intersection_data.w[face_v];
    										
    											for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
    											{
    												b_index = b_start + i_stokes;										

    												// interpolate eta and rho													
													etas[i_stokes] += weight * eta_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index]; 
													rhos[i_stokes] += weight * rho_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index];                                                    

													// interpolate S1 and I1
													S1[i_stokes] += weight * S_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index];	                                                    
													I1[i_stokes] += weight * I_dev.block(i_intersect[face_v] ,j_intersect[face_v],k_intersect[face_v])[b_index];	     												
    											}											
    										}																											                                        

    										K1 = assemble_propagation_matrix_scaled(etas, rhos);                                                                            
    										
    										// optical depth step								
    										dtau = coeff * (eta_I_1 + etas[0]) * intersection_data.distance;									
                                            
    										if (dtau > 0)  std::cout << "ERROR in dtau sign, dtau = " << dtau << std::endl;  
                                            if (dtau == 0) std::cout << "WARNING: dtau = 0, possible e.g. for N_chi = 4"<< std::endl;

    										formal_solver_.one_step(dtau, K1, K2, S1, S2, I1, I2);

                                            // // test
                                            // if (j_theta ==  N_theta - 1 and k_chi == N_chi - 1 and n == 48 and i == i_start and j == j_start)                                                
                                            // {                                                                                                                         
                                            //     // std::cout << "\nk = " << k << std::endl;                                              
                                            //     std::cout << "\nk_global = " << g_dev.global_coord(2, k_aux) << std::endl;                                              
                                            //     // std::cout << "theta = " << theta << std::endl;
                                            //     std::cout << "mu = "  << mu  << std::endl;                                                
                                            //     std::cout << "chi = " << chi << std::endl;                                                
                                            //     // std::cout << "n = "  << n << std::endl;                                                                                            
                                                
                                            //     // std::cout << "dtau = "     << dtau  << std::endl;   
                                            //     // std::cout << "coeff = "    << coeff << std::endl;   
                                            //     // std::cout << "etas[0] = "  << etas[0] << std::endl;  
                                            //     // std::cout << "eta_I_1 = "  << eta_I_1 << std::endl;                                                   
                                            //     // std::cout << "distance = " << intersection_data.distance << std::endl;                                                                                                           

                                            //     // std::cout << "mpi_rank_ = " << mpi_rank_ << std::endl;                                               

                                            //     std::cout << "I1 = "   << I1[0] << std::endl;   
                                            //     // std::cout << "Q1 = "   << I1[1] << std::endl;   
                                            //     // std::cout << "U1 = "   << I1[2] << std::endl;   
                                            //     // std::cout << "V1 = "   << I1[3] << std::endl;   

                                            //     std::cout << "I2 = "   << I2[0] << std::endl;    
                                            //     // std::cout << "Q2 = "   << I2[1] << std::endl;   
                                            //     // std::cout << "U2 = "   << I2[2] << std::endl;                                                   
                                            //     // std::cout << "V2 = "   << I2[3] << std::endl;                                                   
                                                
                                            //     std::cout << "S1[0] = " << S1[0]<< std::endl;
                                            //     std::cout << "S2[0] = " << S2[0]<< std::endl;
                                            //     // std::cout << "S1[3] = " << S1[3]<< std::endl;
                                            //     // std::cout << "S2[3] = " << S2[3]<< std::endl;
                                            //     // std::cout << "S1 = " << std::endl;
                                            //     // for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S1[i_stokes] << std::endl;

                                            //     // std::cout << "S2 = " << std::endl;
                                            //     // for (int i_stokes = 0; i_stokes < 4; ++i_stokes) std::cout << S2[i_stokes] << std::endl;
                                                                            
                                            //     // std::cout << "K1 = " << std::endl;
                                            //     // for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K1[i_stokes] << std::endl;

                                            //     // std::cout << "K2 = " << std::endl;
                                            //     // for (int i_stokes = 0; i_stokes < 16; ++i_stokes) std::cout << K2[i_stokes] << std::endl;
        									// }
                                        
                                            one_step_timer += MPI_Wtime() - start_one;                           									                                        												
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

        }      
          
        if (timing_debug) MPI_Barrier(MPI_COMM_WORLD);
        start_comm = MPI_Wtime();    
        
        I_remap_.from_pblock_to_pgrid(*I_field_serial_, *I_field, tile_number); 

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

    std::list<emission_coefficient_components> components;    

    if (RT_problem_->use_CRD_limit_)
    {
        components.push_back(emission_coefficient_components::epsilon_pCRD_VHP_limit);             
        components.push_back(emission_coefficient_components::epsilon_csc);      

        if (mpi_rank_ == 0) std::cout << "\nUsing CRD emission, components:"<< std::endl;
    }
    else
    {
        // components.push_back(emission_coefficient_components::epsilon_R_II_CONTRIB_FAST);
        components.push_back(emission_coefficient_components::epsilon_R_II_CONTRIB);
        // components.push_back(emission_coefficient_components::epsilon_R_II);
        components.push_back(emission_coefficient_components::epsilon_R_III_GL);
        components.push_back(emission_coefficient_components::epsilon_csc);      

        if (mpi_rank_ == 0) std::cout << "\nUsing PRD emission, components:"<< std::endl;        
    }

    epsilon_fun_ = ecc_sh_ptr_->make_computation_function(components);    

    // Print out emission module
    if (mpi_rank_ == 0) std::cout << ecc_sh_ptr_->emission_components_to_string();    
    
    // module for preconditioner 
    std::list<emission_coefficient_components> components_approx{    
        emission_coefficient_components::epsilon_pCRD_limit,        
        // emission_coefficient_components::epsilon_R_II_AA_PRECOND,
        // emission_coefficient_components::epsilon_R_III,
        // emission_coefficient_components::epsilon_R_II_CONTRIB_EXTREME_FAST,
        emission_coefficient_components::epsilon_csc
    };       
    
    epsilon_fun_approx_ = ecc_sh_ptr_->make_computation_function(components_approx);

    offset_fun_ = rii_include::make_default_offset_function(RT_problem_->N_theta_, RT_problem_->N_chi_, RT_problem_->N_nu_);
        	
	// rii_consts::rii_units::kilometer);		  //

    // set threads number
    // ecc_sh_ptr_->set_threads_number(2);

    // print memory
    unsigned long long b;
    b = ecc_sh_ptr_->bytes();

    if (mpi_rank_ == 0) std::cout << "\n[Memory from set_up_emission_module() = " << (double)b / (1000 * 1024 * 1024) << " GB]" << std::endl;

    auto fs_mem_stat = ecc_sh_ptr_->get_memory_usage_stat();

    if (mpi_rank_ == 0) std::cout << std::endl << fs_mem_stat.to_string() << std::endl;

    if (mpi_rank_ == 0) std::cout << std::endl << fs_mem_stat.sam_memory_stat.to_string() << std::endl;

    //if (mpi_rank_ == 0) std::cout << std::endl << ecc_sh_ptr_->margins_to_string() << std::endl;

    if (mpi_rank_ == 0) std::cout << "--------------------- done -------------------" << std::endl;

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

    PetscInt ix[block_size];

    PetscInt istart, iend; 
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

#ifdef CLOCK_EPSILON
        if (this->mpi_rank_ == 0 and not approx) {
            std::cout << "Start epsilon_computation_function [MAIN], rank: " << this->mpi_rank_ << std::endl;
        }
   
        auto clock = rii_utils::cpu_clock();
        clock.start_clock();
#endif

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

#ifdef CLOCK_EPSILON
	    clock.stop_clock();
        if (this->mpi_rank_ == 0 and not approx) {
              clock.print_clock_h("Execution time of Epsilon emissivity");
        }
#endif


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


// emissivity from scattering (line + continuum)
// emissivity block for a arbitrary direction (mu,chi) is saved in first N_nu block entries of S_field 
// for each spatial point (i,j,k)
void MF_context::update_emission_Omega(const Vec &I_vec, const Real theta, const Real chi)
{      
    const bool include_eps_lth = true;        
    // const bool include_continuum = false;
    const bool include_continuum = RT_problem_->enable_continuum_;

    if (mpi_rank_ == 0)
    {
        std::cout << "\nUpdating emission for theta = " << theta << " and chi = " << chi << std::endl;

        if (include_continuum)
        {
            std::cout << "Including continuum" << std::endl;
        }
        else
        {
            std::cout << "NOT including continuum" << std::endl;
        }
    }

    PetscErrorCode ierr; 
           
    const auto g_dev       = RT_problem_->space_grid_      ->view_device();  
    const auto eta_dev     = RT_problem_->eta_field_Omega_ ->view_device();
    const auto S_dev_Omega = RT_problem_->S_field_Omega_   ->view_device(); 
    
    const auto block_size = RT_problem_->block_size_;   
    const auto N_nu       = RT_problem_->N_nu_;     
    
    // field range indeces 
    const int i_start = g_dev.margin[0]; 
    const int j_start = g_dev.margin[1];
    const int k_start = g_dev.margin[2];

    const int i_end = i_start + g_dev.dim[0]; 
    const int j_end = j_start + g_dev.dim[1];
    const int k_end = k_start + g_dev.dim[2];   

    std::vector<double>  input(block_size);        

    PetscInt ix[block_size];

    PetscInt istart, iend; 
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
        const std::string scattering_model = (RT_problem_->use_CRD_limit_) ? "CRD" : "PRD";        
        auto epsilon_computation_Omega = ecc_sh_ptr_->make_computation_function_arbitrary_direction(scattering_model, include_eps_lth, include_continuum);

        ecc_sh_ptr_->update_incoming_field(i, j, k, offset_fun_, input.data());

        // get IQUV for (theta, chi direction)
        auto IQUV_matrix_sh_ptr = epsilon_computation_Omega(i, j, k, theta, chi);           
                
        // update S_field_ from output scaling by eta_I
        double eta_I_inv; 

        // index
        int b;        

        for (int n_nu = 0; n_nu < N_nu; n_nu++)
        {
            b = 4 * n_nu;

            eta_I_inv = 1.0 / (eta_dev.block(i,j,k)[b]);            

            for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
            {
                S_dev_Omega.block(i,j,k)[b + i_stokes] = eta_I_inv * (*IQUV_matrix_sh_ptr)(n_nu, i_stokes);                
            }                                    
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



void MF_context::init_serial_fields(const int n_tiles){
    
    auto block_size = RT_problem_->block_size_;

    auto N_x = RT_problem_->N_x_;
    auto N_y = RT_problem_->N_y_;
    auto N_z = RT_problem_->N_z_;

    auto N_theta = RT_problem_->N_theta_;
    auto N_chi   = RT_problem_->N_chi_;
    auto N_nu    = RT_problem_->N_nu_;

    // set the number of local rays and tiles
    n_tiles_ = n_tiles;    
    n_local_rays_ = block_size/mpi_size_;

    if (n_local_rays_ < 4) 
    {
        if (mpi_rank_ == 0) std::cout << "WARNING: mpi_size > number of rays" << std::endl;
        n_local_rays_ = 4; 
    } 
    else
    {
        if (n_local_rays_ * mpi_size_ != block_size){ 
            std::cout << "ERROR file: " << __FILE__ << " line: " << __LINE__ << std::endl;
            std::cout << "ERROR in init_serial_fields(): block_size/mpi_size_ not integer" << std::endl;
            std::cout << "ERROR block_size = " << block_size << std::endl;
            std::cout << "ERROR mpi_size = " << mpi_size_ << std::endl;
            std::cout << "ERROR n_local_rays_ = " << n_local_rays_ << std::endl;
            std::cout << "ERROR block_size % mpi_size_ = " << (block_size % mpi_size_) << std::endl;
        }  
    }
    
    const int N_thea_chi = N_theta * N_chi;

    // check block decomposition
    if (mpi_size_ < N_thea_chi * N_nu)
    {
        if (mpi_size_ > N_thea_chi)
        {        
            if (( N_thea_chi * N_nu / mpi_size_) * mpi_size_ != N_thea_chi * N_nu)
            {
                std::stringstream ss;
                ss << "ERROR: with block decomposition I, at line  " << __LINE__ << " in file " <<  __FILE__ << std::endl;
                ss << "ERROR: N_thea_chi = " << N_thea_chi << ", N_nu = " << N_nu << ", mpi_size_ = " << mpi_size_;
                ss << ", N_thea_chi * N_nu = " << N_thea_chi * N_nu << ", ( N_thea_chi * N_nu / mpi_size_) * mpi_size_ = " << ( N_thea_chi * N_nu / mpi_size_) * mpi_size_;
                throw std::runtime_error(ss.str());
            } 
        }
        else if (mpi_size_ > N_theta)
        {
            if (( N_thea_chi / mpi_size_) * mpi_size_ != N_thea_chi)
            {       
                std::stringstream ss;
                ss << "ERROR: with block decomposition II, at line  " << __LINE__ << " in file " <<  __FILE__ << std::endl;
                ss << "ERROR: N_thea_chi = " << N_thea_chi << ", mpi_size_ = " << mpi_size_;
                ss << ", N_thea_chi / mpi_size_ = " << N_thea_chi / mpi_size_ << ", ( N_thea_chi / mpi_size_) * mpi_size_ = " << ( N_thea_chi / mpi_size_) * mpi_size_;     
                throw std::runtime_error(ss.str());
            } 
        }
    }

    if (n_local_rays_ % 4 != 0) std::cout << "ERROR in init_serial_fields(): n_local_rays_ should be divisible by 4" << std::endl;        

    tile_size_ = n_local_rays_/n_tiles_;

    if (tile_size_ * n_tiles_ != n_local_rays_) std::cout << "ERROR in init_serial_fields(): n_local_rays_/n_tiles_ not integer" << std::endl;        
    if (tile_size_ % 4 != 0)                    std::cout << "ERROR in init_serial_fields(): tile_size_ should be divisible by 4" << std::endl;            

    // init serial grid
    const bool use_ghost_layers = false;

    space_grid_serial_ = std::make_shared<Grid_t>();    
    space_grid_serial_->init(MPI_COMM_SELF, {N_x, N_y, N_z}, {1, 1, 0}, {}, use_ghost_layers); 

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


void MF_context::init_serial_fields_Omega(){
    
    auto N_z  = RT_problem_->N_z_;
    auto N_nu = RT_problem_->N_nu_;
    auto block_size = 4 * N_nu;

    // set the number of local rays and tiles   
    local_block_size_ = block_size/mpi_size_;        

    if (local_block_size_ < 4)
    {
        if (mpi_rank_ == 0) std::cout << "WARNING: mpi_size > number of rays" << std::endl;
        
        local_block_size_ = 4; 
        
        // // some ranks can stay idle
        // if (mpi_rank_ >= N_nu)
        // {
        //     idle_processor_Omega_ = true;         TODO  local_block_size_ = 1 here?
        // }
    } 
    else
    {
        if (local_block_size_ * mpi_size_ != block_size) { 
            std::cout << "ERROR file: " << __FILE__ << " line: " << __LINE__ << std::endl;
            std::cout << "ERROR in init_serial_fields(): block_size/mpi_size_ not integer" << std::endl;
            std::cout << "ERROR block_size = " << block_size << std::endl;
            std::cout << "ERROR mpi_size = " << mpi_size_ << std::endl;
            std::cout << "ERROR n_local_rays_ = " << n_local_rays_ << std::endl;
            std::cout << "ERROR block_size % mpi_size_ = " << (block_size % mpi_size_) << std::endl;
        }    
    }
    
    if (local_block_size_ % 4 != 0) std::cout << "ERROR in init_serial_fields(): local_block_size_ should be divisible by 4" << std::endl;        
        
    // create serial fields 
    I_field_serial_Omega_   = std::make_shared<Field_t>("I_serial", space_grid_serial_, local_block_size_); 
    S_field_serial_Omega_   = std::make_shared<Field_t>("S_serial", space_grid_serial_, local_block_size_);    

    eta_field_serial_Omega_ = std::make_shared<Field_t>("eta_serial", space_grid_serial_, local_block_size_); // here could tiles also be used to reduce mem footprint
    rho_field_serial_Omega_ = std::make_shared<Field_t>("rho_serial", space_grid_serial_, local_block_size_);

    // allocate
    I_field_serial_Omega_  ->allocate_on_device();     
    S_field_serial_Omega_  ->allocate_on_device();     

    eta_field_serial_Omega_->allocate_on_device();     
    rho_field_serial_Omega_->allocate_on_device();     

    // init remaps 
    I_remap_Omega_.init(*(RT_problem_->I_field_Omega_), *I_field_serial_Omega_); 
    S_remap_Omega_.init(*(RT_problem_->S_field_Omega_), *S_field_serial_Omega_);      

    // apply BC on I_Field
    const auto I_Omega_dev = RT_problem_->I_field_Omega_->view_device();      
    const auto W_T_dev     = RT_problem_->W_T_          ->view_device();     
    const auto g_dev       = RT_problem_->space_grid_   ->view_device();  

    sgrid::parallel_for("APPLY BC", RT_problem_->space_grid_->md_range(), SGRID_LAMBDA(int i, int j, int k) {
                                    
        // just in max depth
        if (g_dev.global_coord(2, k) == (N_z - 1))        
        {       
            const Real W_T_deep = W_T_dev.ref(i,j,k);            
            
            for (int b = 0; b < block_size; b = b + 4) 
            {
                I_Omega_dev.block(i,j,k)[b] = W_T_deep;                                
            }                       
        }
    });     

    // init BC in serial grid
    I_remap_Omega_.from_pgrid_to_pblock(*RT_problem_->I_field_Omega_, *I_field_serial_Omega_, 0);                                    
} 


void RT_solver::print_info(){

    // print some output
    if (mpi_rank_ == 0)
    {        
        if (mf_ctx_.use_single_long_step_)
        {
            std::cout << "Using a single step for long rays." << std::endl;
        }
        else
        {
            std::cout << "Using multiple steps for long rays." << std::endl;
        }

        if (mf_ctx_.use_always_long_ray_) std::cout << "Only long rays are used." << std::endl;

        std::cout << "\n========= Serial formal solver parameters =========" << std::endl;
        std::cout << "n_local_rays = " << mf_ctx_.n_local_rays_ << " (block_size/mpi_size)" << std::endl;             
        std::cout << "tile_size    = " << mf_ctx_.tile_size_    << " (n_local_rays/n_tiles)" << std::endl;    
        std::cout << "n_tiles      = " << mf_ctx_.n_tiles_ << std::endl;                          
        std::cout << "===================================================" << std::endl;
    } 
}


void RT_solver::assemble_rhs(){ 

    // with test = true data structures are created but not filled
    const bool test = false;

  	if (mpi_rank_ == 0 and (not test)) std::cout << "\n++++++ Assembling right hand side...+++++++++";
    if (mpi_rank_ == 0 and test)       std::cout << "\n+++++++++++ RHS TEST RHS TEST +++++++++++++";
 
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

    		std::vector<int> local_idx;

    		for (int b = 0; b < block_size; b++) 
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
    mf_ctx_->formal_solve_global(RT_problem->I_field_, RT_problem->S_field_, 0.0); 
      
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
   
    if (RT_problem->mpi_rank_ == 0) printf("Update preconditioner emission:\t\t%g seconds\n", MPI_Wtime() - start);              

    // fill rhs_ from formal solve with zero bc     
    mf_ctx_->formal_solve_global(RT_problem->I_field_, RT_problem->S_field_, 0.0);
    
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
    PetscInt iterations;
    ierr = KSPGetIterationNumber(mf_ctx->pc_solver_, &iterations);CHKERRQ(ierr);
    if (mf_ctx->mpi_rank_ == 0) std::cout << "Preconditioner iterations: " << iterations << std::endl;

	return ierr;
}
