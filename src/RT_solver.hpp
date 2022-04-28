#ifndef RT_solver_hpp
#define RT_solver_hpp

#include "Formal_solver.hpp"
#include "RT_problem.hpp"
#include <rii_emission_coefficient_3D.h>
#include "sgrid_ReMap.hpp"

extern PetscErrorCode UserMult(Mat mat,Vec x,Vec y);
extern PetscErrorCode UserMult_approx(Mat mat,Vec x,Vec y);
extern PetscErrorCode MF_pc_Destroy(PC pc);
extern PetscErrorCode MF_pc_Apply(PC pc,Vec x,Vec y);


// struct for ray - grid intersection 
typedef struct t_intersect {
    int ix[4], iy[4], iz[4];
    double w[4];
    double distance;
} t_intersect;


typedef struct t_xyinters {
    int plane; // 0 = const x, 1 = const y, 2 = const z
    int ix, iy, iz;
    double x, y, z;
} t_xyinters;


// matrix-free (MF) structure
struct MF_context {

	std::shared_ptr<RT_problem> RT_problem_;	

	Formal_solver formal_solver_;

	// preconditioner data structures 
	KSP pc_solver_;

	// MPI varables
	int mpi_rank_;
	int mpi_size_;
	
	// if false linear interpolation is used
	bool use_log_interpolation_ = false;
	bool use_single_long_step_  = true;
	bool use_always_long_ray_   = true;

	// serial objects for formal solution
	Grid_ptr_t  space_grid_serial_;	

	sgrid::ReMap<Field_t> I_remap_;
	sgrid::ReMap<Field_t> S_remap_;
	
	Field_ptr_t I_field_serial_;
	Field_ptr_t S_field_serial_;
	Field_ptr_t eta_field_serial_;
	Field_ptr_t rho_field_serial_;

	// total number of rays a single processor will handle, n_local_rays_ = block_size/mpi_size
	int n_local_rays_;

	// number of tiles each processor will handle, n_tiles_ = 1 default, to be increased to reduce memmory usage
    int n_tiles_;   
		
	// number of rays a single processor will handle at one time, tile_size_ = n_local_rays_/n_tiles_
	int tile_size_;	
	
	// pointer for emission module and offset
	std::shared_ptr<rii_include::emission_coefficient_computation_3D> ecc_sh_ptr_;
	rii_include::emission_coefficient_computation_3D::compute_node_3D_function_type epsilon_fun_; 
	rii_include::emission_coefficient_computation_3D::compute_node_3D_function_type epsilon_fun_approx_; 
	rii_include::offset_function_cartesian offset_fun_;	
	
	// change data format
	void field_to_vec(const Field_ptr_t field, Vec &v);
	void vec_to_field(Field_ptr_t field, const Vec &v);
	
	// find intersection
	void find_intersection(double theta, double chi, const double Z_down, const double Z_top, const double L, t_intersect *T);
	std::vector<t_intersect> find_prolongation(double theta, double chi, const double dz, const double L);
	std::vector<double> long_ray_steps(const std::vector<t_intersect> T, const Field_ptr_t I_field, const Field_ptr_t S_field, const int i, const int j, const int k, const int block_index);
	std::vector<double> single_long_ray_step(const std::vector<t_intersect> T, const Field_ptr_t I_field, const Field_ptr_t S_field, const int i, const int j, const int k, const int block_index);

	void get_2D_weigths(const double x, const double y, double *w);

	// formal solver	
	void formal_solve_global(Field_ptr_t I_field, Field_ptr_t I_field_serial, const Field_ptr_t S_field, const Field_ptr_t S_field_serial, const Real I0);		
	
	void apply_bc(Field_ptr_t I_field, const Real I0);	
		
	// emission module from Simone
	void set_up_emission_module();
		
	// update emission in all spatial points (given the current I_field_, update S_field_)
	void update_emission(const Vec &I_field, const bool approx = false);

	// init serial fields (serial eta and rho are filled)
	void init_serial_fields(const size_t n_tiles);	
};

class RT_solver
{
public:
	RT_solver(const std::shared_ptr<RT_problem> RT_problem, input_string formal_solver = "implicit_Euler", const bool using_prec = true) 
	{		
		PetscErrorCode ierr;

		// assign MPI varaibles and init mf_ctx_
    	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
    	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);  

    	RT_problem_ = RT_problem;  
    	using_prec_ = using_prec;    	

    	mf_ctx_.RT_problem_    = RT_problem;  
    	mf_ctx_.mpi_rank_      = mpi_rank_;
    	mf_ctx_.mpi_size_      = mpi_size_;   
    	mf_ctx_.formal_solver_ = Formal_solver(formal_solver);    

    	// init serial grids for formal solution
    	const size_t n_tiles = 1; // TODO: now fixed
    	mf_ctx_.init_serial_fields(n_tiles);	
    	    	
    	mf_ctx_.set_up_emission_module();  	  

    	// print some output 
    	print_info();

    	// assemble rhs
    	assemble_rhs();
    	// save_vec(rhs_, "../output/rhs.m" ,"rhs_3d");          	
  
    	// set linear system
		int local_size = RT_problem_->local_size_;
		ierr = VecGetLocalSize(rhs_, &local_size);CHKERRV(ierr); 		

		// init user defined Mat mult
		ierr = MatCreateShell(PETSC_COMM_WORLD,local_size,local_size,RT_problem_->tot_size_,RT_problem_->tot_size_,(void*)&mf_ctx_,&MF_operator_);CHKERRV(ierr); 
		ierr = MatShellSetOperation(MF_operator_,MATOP_MULT,(void(*)(void))UserMult);CHKERRV(ierr);

    	// set Krylov solver
    	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp_solver_);CHKERRV(ierr);
    	ierr = KSPSetOperators(ksp_solver_,MF_operator_,MF_operator_);CHKERRV(ierr);	    		
    	ierr = KSPSetType(ksp_solver_,ksp_type_);CHKERRV(ierr); 

    	// set preconditioner
    	ierr = KSPGetPC(ksp_solver_,&pc_);CHKERRV(ierr);    		    	    		

    	if (using_prec_)
    	{    	    		
			// set MF_operator_approx_		
			ierr = MatCreateShell(PETSC_COMM_WORLD,local_size,local_size,RT_problem_->tot_size_,RT_problem_->tot_size_,(void*)&mf_ctx_,&MF_operator_approx_);CHKERRV(ierr); 
			ierr = MatShellSetOperation(MF_operator_approx_,MATOP_MULT,(void(*)(void))UserMult_approx);CHKERRV(ierr);		

			// set PC solver 
			ierr = KSPCreate(PETSC_COMM_WORLD,&mf_ctx_.pc_solver_);CHKERRV(ierr);
    		ierr = KSPSetOperators(mf_ctx_.pc_solver_,MF_operator_approx_,MF_operator_approx_);CHKERRV(ierr);	    		
    		ierr = KSPSetType(mf_ctx_.pc_solver_,KSPGMRES);CHKERRV(ierr); 

    		// const int max_its = 10;
    		// // const double r_tol = 1e-3;
    		// // // ierr = KSPSetFromOptions(mf_ctx_.pc_solver_);CHKERRV(ierr);    		
    		// // ierr = KSPSetTolerances(mf_ctx_.pc_solver_,r_tol,PETSC_DEFAULT,PETSC_DEFAULT, PETSC_DEFAULT);CHKERRV(ierr);
    		// ierr = KSPSetTolerances(mf_ctx_.pc_solver_,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT, max_its);CHKERRV(ierr);
    		// ierr = KSPSetNormType(mf_ctx_.pc_solver_, KSP_NORM_NONE);CHKERRV(ierr);


    		// set PC
    		ierr = PCSetType(pc_,PCSHELL);CHKERRV(ierr);
			ierr = PCShellSetContext(pc_, &mf_ctx_);CHKERRV(ierr);		
			ierr = PCShellSetApply(pc_,MF_pc_Apply);CHKERRV(ierr);				
			ierr = PCShellSetDestroy(pc_,MF_pc_Destroy);CHKERRV(ierr);	
    	}
    	else
    	{
    		ierr = PCSetType(pc_,PCNONE);CHKERRV(ierr);
    	}

    	// extra options from command line   	
    	ierr = KSPSetFromOptions(ksp_solver_);CHKERRV(ierr);
    	// ierr = PCSetFromOptions(pc_);CHKERRV(ierr);	     

		// test
		// ierr = MatMult(MF_operator_, rhs_, RT_problem_->I_vec_);CHKERRV(ierr);	        				
		// const std::string filename =  "../output/rhs_" + std::to_string(mpi_size_) + ".m";
  //   	const std::string varible  =  "rhs" + std::to_string(mpi_size_);
  //   	save_vec(rhs_, filename.c_str(), varible.c_str());                 	  		
	}

	// solve linear system
	inline void solve()
	{	
		PetscErrorCode ierr;

		Real start = MPI_Wtime();	
				
		if (mpi_rank_ == 0) std::cout << "\nStart linear solve..." << std::endl;			
		ierr = KSPSolve(ksp_solver_, rhs_, RT_problem_->I_vec_);CHKERRV(ierr);

		MPI_Barrier(MPI_COMM_WORLD); Real end = MPI_Wtime();
		if (mpi_rank_ == 0) std::cout << "Solve time (s) = " << end - start << std::endl;	

		// update I_field for later use
		mf_ctx_.vec_to_field(RT_problem_->I_field_, RT_problem_->I_vec_);
	}

	inline void apply_formal_solver()
	{		
		Real start = MPI_Wtime();		

		// // set source fun
		// auto S_dev = RT_problem_->S_field_->view_device();

	 //    sgrid::parallel_for("INIT S", RT_problem_->space_grid_->md_range(), KOKKOS_LAMBDA(int i, int j, int k) 
	 //    {         
	 //        auto *block = S_dev.block(i, j, k);
	         
	 //        for (int b = 0; b < (int)RT_problem_->block_size_; ++b) 
	 //        {
	 //        	block[b] = 1e-20;        	
	 //        }
	 //    });

		if (mpi_rank_ == 0) std::cout << "Start formal solve..." << std::endl;

		mf_ctx_.formal_solve_global(RT_problem_->I_field_, mf_ctx_.I_field_serial_, RT_problem_->S_field_, mf_ctx_.S_field_serial_, 1.0);		
						
		MPI_Barrier(MPI_COMM_WORLD); Real end = MPI_Wtime();
		if (mpi_rank_ == 0) std::cout << "Formal solve time (s) = " << end - start << std::endl;	

		// // update I_vec for later use
		// mf_ctx_.field_to_vec(RT_problem_->I_field_, RT_problem_->I_vec_);
	}	

	inline void compute_emission()
	{
		Real start = MPI_Wtime();		

		if (mpi_rank_ == 0) std::cout << "Computing emission..." << std::endl;
		
		// compute new emission in S_field_ 
  		mf_ctx_.update_emission(RT_problem_->I_vec_);   

    	// exchange updated S_field_             
    	RT_problem_->S_field_->exchange_halos();   	
			
		MPI_Barrier(MPI_COMM_WORLD); Real end = MPI_Wtime();
		if (mpi_rank_ == 0) std::cout << "Computing emission took (s) = " << end - start << std::endl;	

		mf_ctx_.field_to_vec(RT_problem_->S_field_, RT_problem_->I_vec_);
	}
	

	inline void test_transfer()
	{			
		PetscErrorCode ierr; 

		const auto g_dev      = RT_problem_->space_grid_->view_device();
		const auto field_dev  = RT_problem_->I_field_->view_device();	
		const auto block_size = RT_problem_->block_size_;

		// indeces
		const int i_start = g_dev.margin[0]; 
		const int j_start = g_dev.margin[1];
		const int k_start = g_dev.margin[2];

		const int i_end = i_start + g_dev.dim[0];
		const int j_end = j_start + g_dev.dim[1];
		const int k_end = k_start + g_dev.dim[2];

		int counter = 0;

		int istart;	
		ierr = VecGetOwnershipRange(RT_problem_->I_vec_, &istart, NULL);CHKERRV(ierr);	

		// init
		for (int k = k_start; k < k_end; ++k)					
		{															
			for (int j = j_start; j < j_end; ++j)
			{
				for (int i = i_start; i < i_end; ++i)				
				{
					for (int b = 0; b < (int)block_size; b++) 
					{			
						field_dev.block(i, j, k)[b] = istart + counter;							

						counter++;
					}							
				}
			}
		}

		// test 1
		mf_ctx_.field_to_vec(RT_problem_->I_field_, RT_problem_->I_vec_);

		save_vec(RT_problem_->I_vec_,  "../output/vec1.m" ,"vec_1");   

		// test 2 
		mf_ctx_.vec_to_field(RT_problem_->I_field_, RT_problem_->I_vec_);
		mf_ctx_.field_to_vec(RT_problem_->I_field_, RT_problem_->I_vec_);

		save_vec(RT_problem_->I_vec_,  "../output/vec2.m" ,"vec_2");   

		// const std::string filename =  "../output/I_" + std::to_string(mpi_size_) + ".m";
	}

private:	

	// MPI varables
	int mpi_rank_;
	int mpi_size_;

	// for astmospheric data and parallel fields
	std::shared_ptr<RT_problem> RT_problem_;	
	
	// MF context
	MF_context mf_ctx_;
	
	// linear system quantities
	Mat MF_operator_;
	Mat MF_operator_approx_;
	Vec rhs_;
	
	KSP ksp_solver_;
	KSPType ksp_type_ = KSPGMRES;
	PC pc_;
	
	bool using_prec_;
			
	// assemble Lam[eps_th] + t
	void assemble_rhs();		

	void print_info();	
};


#endif 
// 
