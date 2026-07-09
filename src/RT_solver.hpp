#ifndef RT_solver_hpp
#define RT_solver_hpp

#include "Formal_solver.hpp"
#include "RT_problem.hpp"
#include "GridManager/GridManager.hpp"

extern PetscErrorCode UserMult(Mat mat,Vec x,Vec y);
extern PetscErrorCode UserMult_approx(Mat mat,Vec x,Vec y);
extern PetscErrorCode UserMult_JKQ(Mat mat,Vec x,Vec y);
extern PetscErrorCode MF_pc_Destroy(PC pc);
extern PetscErrorCode MF_pc_Apply(PC pc,Vec x,Vec y);

unsigned int get_RII_contrib_block_size();

void set_RII_contrib_block_size(const unsigned int block_size);

// structs for ray - grid intersection 
typedef struct t_intersect {
    int ix[4], iy[4], iz[4];
    double w[4];
    double distance;
} t_intersect;

enum t_plane {
    I_YZ = 0,
    I_XZ,
    I_XY
};

typedef struct t_xyzinters {
    t_plane plane;
    double x, y, z; // point of intersection
    int ix, iy, iz; // left-back-from index of the intersected segment
    double t; // distance from the origin
} t_xyzinters;


// matrix-free (MF) structure
struct MF_context {

	std::shared_ptr<RT_problem> RT_problem_;

	Formal_solver formal_solver_;
	Formal_solver formal_solver_unpol_;

	// for preconditioning
	bool approx_formal_solver_; // TODO put this into config

	// preconditioner data structures 
	KSP pc_solver_;

	// MPI varables
	int mpi_rank_;
	int mpi_size_;
	
	const bool use_long_characteristics_ = true; // gor rays intersecting vertical planes
	const bool use_single_long_step_     = false; 

	// use JKQ vector as unknown
	bool ksp_use_J_KQ_;
	int tot_J_KQ_size_; // e.g. N_s * J_KQ_size_
	int local_J_KQ_size_; 
	const int J_KQ_size_ = 18;

	// reduced models flags for preconditioner
	const bool unpolarized_prec_ = false; // TEST
	bool pc_use_J_KQ_; 
	
	// formal solution in arbitrary direction
	bool formal_solution_Omega_ = false;
	
	// serial objects for formal solution
	Grid_ptr_t space_grid_serial_;	

	ReMap3D Pol_remap;
	
	Field_ptr_t I_field_serial_;
	Field_ptr_t S_field_serial_;
	Field_ptr_t eta_field_serial_;
	Field_ptr_t rho_field_serial_;

	std::vector<Real> W_T_ij_serial_; 

	// auxiliary unpolarized vecotors
	Field_ptr_t   I_unpol_field_serial_;
	Field_ptr_t   S_unpol_field_serial_;

	Vec x_unpol_, y_unpol_, x_pol_;
	IS pol_indeces_ = nullptr;

	// auxiliary J_KQ vecotors
	Vec x_J_KQ_, y_J_KQ_;

	ReMap3D Unpol_remap;

	// data structures a single direction Omega (if needed)
	ReMap3D Omega_remap;
	
	Field_ptr_t I_field_serial_Omega_;
	Field_ptr_t S_field_serial_Omega_;
	Field_ptr_t eta_field_serial_Omega_;
	Field_ptr_t rho_field_serial_Omega_;

	// total number of rays a single processor will handle, n_local_rays_ = block_size/mpi_size
	int n_local_rays_;
	int n_local_rays_unpol_;
	int local_block_size_;

	// number of tiles each processor will handle, n_tiles_ = 1 default, to be increased to reduce memmory usage
    int n_tiles_;   
		
	// number of rays a single processor will handle at one time, tile_size_ = n_local_rays_/n_tiles_
	int tile_size_;	
	
	// pointer for emission module and offset
	std::shared_ptr<rii_include::emission_coefficient_computation_3D> ecc_sh_ptr_;
	rii_include::emission_coefficient_computation_3D::compute_node_3D_function_type epsilon_fun_; 
	rii_include::emission_coefficient_computation_3D::compute_node_3D_function_type epsilon_fun_approx_; 
	rii_include::emission_coefficient_computation_3D::compute_node_3D_function_type epsilon_fun_csc_; 
	rii_include::emission_coefficient_computation_3D::compute_node_3D_function_eps_from_JKQ epsilon_fun_J_KQ_; 
	rii_include::emission_coefficient_computation_3D::compute_node_3D_function_JKQ_vals_type compute_JKQ_values_;
	rii_include::offset_function_cartesian offset_fun_;	

#if ACC_SOLAR_3D == _ON_
	rii_include::emission_coefficient_computation_3D::start_device_handler_function_type start_device_handler_fun_;
#endif
	
	// change data format
	void field_to_vec(const Field_ptr_t field, Vec &v, const int block_size = -1); //UNUSED
	void vec_to_field(Field_ptr_t field, const Vec &v, const int block_size = -1); //UNUSED

	void I_vec_to_J_KQ_vec(const Vec &I_vec, Vec &J_KQ_vec);
	void I_field_to_J_KQ_vec(const Field_ptr_t field, Vec &J_KQ_vec);
	
	// find intersection
	void find_intersection(double theta, double chi, const double Z_down, const double Z_top, const double L, t_intersect *T);
	std::vector<t_intersect> find_prolongation(double theta, double chi, const double dz, const double L);
	std::vector<double> long_ray_steps(const std::vector<t_intersect> T, const Field_ptr_t I_field, const Field_ptr_t S_field, const int i, const int j, const int k, const int block_index);
	std::vector<double> long_ray_steps_quadratic(const std::vector<t_intersect> T, const Field_ptr_t I_field, const Field_ptr_t S_field, const int i, const int j, const int k, const int block_index, bool print_flag);
	std::vector<double> single_long_ray_step(const std::vector<t_intersect> T, const Field_ptr_t I_field, const Field_ptr_t S_field, const int i, const int j, const int k, const int block_index);

	// for the unpolarized case
	double long_ray_steps_unpolarized(          const std::vector<t_intersect> T, const Field_ptr_t I_field, const Field_ptr_t S_field, const int i, const int j, const int k, const int block_index, bool print_flag);
	double long_ray_steps_quadratic_unpolarized(const std::vector<t_intersect> T, const Field_ptr_t I_field, const Field_ptr_t S_field, const int i, const int j, const int k, const int block_index, bool print_flag);

	void get_2D_weigths(const double x, const double y, double *w);

	// formal solver	
	void formal_solve_global(Field_ptr_t I_field, const Field_ptr_t S_field, const Real I0, const bool approx_formal_solver = false);		
	void formal_solve_ray(const double mu, const double chi);		
	void formal_solve_unpolarized(Field_ptr_t I_field, const Field_ptr_t S_field, const Real I0);		
	void formal_solve_1_5D(Field_ptr_t I_field, const Field_ptr_t S_field, const Real I0);	
	void formal_solve(Field_ptr_t I_field, const Field_ptr_t S_field, const Real I0, const bool approx_formal_solver = false);		

	// formation height of a specific LOS
	void get_formation_height(const double mu, const double chi);

	void apply_bc(       Field_ptr_t I_field, const Real I0, const bool polarized = true);	
	void apply_bc_serial(Field_ptr_t I_field, const Real I0, const bool polarized = true);	
		
	// emission module from Simone
	void set_up_emission_module();
		
	// update emission in all spatial points (given the current I_field_, update S_field_)
	void update_emission(const Vec &I_field, const bool approx = false);
	void update_emission_J_KQ(const Vec &I_field);

	// update emission in all spatial points (given the current I_field_, update S_field_) for an arbitrary direction 
	void update_emission_Omega(const Vec &I_field, const double theta, const double chi);

	// init serial fields (serial eta and rho are filled)
	void init_serial_fields(const int n_tiles);	
	void init_unpol_fields();
	void init_J_KQ_vectors();

	// init serial fields (serial eta and rho are filled) in a single direction
	void init_serial_fields_Omega();		

	// moving from polarized and unpolarized vetors
	void polarized_to_unpolarized(Vec &pol_v, Vec &unpol_v);
	void unpolarized_to_polarized(Vec &unpol_v, Vec &pol_v);
};

static PetscErrorCode precond_ksp_monitor(KSP ksp, PetscInt it, PetscReal rnorm, void *ctx)
{
	(void)ctx;
	Vec r;
	PetscReal truenorm;
	PetscCall(KSPBuildResidual(ksp, NULL, NULL, &r));
	PetscCall(VecNorm(r, NORM_2, &truenorm));
	PetscCall(VecDestroy(&r));
	return PetscPrintf(PETSC_COMM_WORLD, "%3d PRECOND Residual norm %14.12e  ||pre_r(i)||/||pre_b|| %14.12e\n", (int)it, (double)rnorm, (double)truenorm);
}

class RT_solver
{
public:
	RT_solver(const std::shared_ptr<RT_problem> RT_problem) 
	{		
		PetscErrorCode ierr;

		// assign MPI varaibles and init mf_ctx_
    	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
    	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);  

    	RT_problem_ = RT_problem;  
    	auto cfg    = RT_problem_->cfg_;

    	// load input options 
    	using_prec_  = cfg.use_prec;    	    	
    	ksp_type_    = cfg.solver.ksp_solver_type;
    	pc_ksp_type_ = cfg.prec.pc_solver_type;	    	

    	mf_ctx_.RT_problem_ = RT_problem;  
    	mf_ctx_.mpi_rank_   = mpi_rank_;
    	mf_ctx_.mpi_size_   = mpi_size_;   

    	mf_ctx_.formal_solver_ = Formal_solver(cfg.formal_solver);         	

    	mf_ctx_.approx_formal_solver_ = cfg.prec.pc_formal_solver_approx;

    	mf_ctx_.ksp_use_J_KQ_ = cfg.solver.ksp_use_J_KQ;
    	mf_ctx_.pc_use_J_KQ_  = cfg.prec.pc_use_J_KQ;

    	// safety check when J_KQ solution is enabled
    	const std::string emissivity_model = emissivity_model_to_string_long(cfg.emissivity_model);
    	
		if (mf_ctx_.ksp_use_J_KQ_ and emissivity_model.find("PRD") != std::string::npos)
		{
		    if (mpi_rank_ == 0)
		    {
		        std::cout << "WARNING: J_KQ solver cannot be used with PRD "
		                  << "(emissivity_model = " << emissivity_model << "). "
		                  << "Switching to ksp_use_J_KQ_ = false.\n";
		    }
		    
		    mf_ctx_.ksp_use_J_KQ_ = false;
		}   	
    
		if (cfg.formal_solver == "BESSER" and emissivity_model.find("PRD") != std::string::npos and strcmp(ksp_type_, KSPRICHARDSON) != 0 and mpi_rank_ == 0)
                {
                   std::cout << "WARNING: BESSER with GMRES can stagnate in PRD, switch to DELO_linear or KSPRICHARDSON for better stability.\n";
                }

    	// init serial grids for formal solution (not needed for 1.5D)
    	if (RT_problem_->use_1_5D_approx_)
    	{
    		mf_ctx_.n_local_rays_ = RT_problem_->block_size_;
    	}
    	else
    	{
    		const int n_tiles = 1; 
	    	mf_ctx_.init_serial_fields(n_tiles);	
    	}
    	    	
    	// init unpolarized formal solver and data structures
    	if (mf_ctx_.pc_use_J_KQ_ or mf_ctx_.ksp_use_J_KQ_)    	
    	{    		
    		mf_ctx_.init_J_KQ_vectors();    		
    	}
    	else if (mf_ctx_.unpolarized_prec_) 
    	{
    		mf_ctx_.RT_problem_->allocate_unpolarized_fields();
    		mf_ctx_.init_unpol_fields();    		
    		mf_ctx_.formal_solver_unpol_ = Formal_solver("SC_parabolic"); // or SC_parabolic 		
    	}
    	   
    	mf_ctx_.set_up_emission_module();  	      	    	
        	
    	// assemble rhs
    	assemble_rhs();
    	// save_vec(rhs_, "../output/rhs.m" ,"rhs_3d");  
    	
    	// //test
    	// mf_ctx_.field_to_vec(RT_problem_->eta_field_, RT_problem_->I_vec_);     
    	// save_vec(RT_problem_->I_vec_, "../output/eta_t.m" ,"etat");     
    	
    	// print some output 
    	if (RT_problem_->verbose_) print_info();
  
    	// set linear system
		PetscInt local_size = RT_problem_->local_size_;
		ierr = VecGetLocalSize(rhs_, &local_size);CHKERRV(ierr); 		

		// init user defined Mat mult
		if (mf_ctx_.ksp_use_J_KQ_)
		{
			ierr = MatCreateShell(PETSC_COMM_WORLD,mf_ctx_.local_J_KQ_size_,mf_ctx_.local_J_KQ_size_,mf_ctx_.tot_J_KQ_size_,mf_ctx_.tot_J_KQ_size_,
													   (void*)&mf_ctx_,&MF_operator_);CHKERRV(ierr); 	
			ierr = MatShellSetOperation(MF_operator_,MATOP_MULT,(void(*)(void))UserMult_JKQ);CHKERRV(ierr);
		}
		else // standard solve with [I Q U V] unknown
		{			
			ierr = MatCreateShell(PETSC_COMM_WORLD,local_size,local_size,RT_problem_->tot_size_,RT_problem_->tot_size_,(void*)&mf_ctx_,&MF_operator_);CHKERRV(ierr); 
			ierr = MatShellSetOperation(MF_operator_,MATOP_MULT,(void(*)(void))UserMult);CHKERRV(ierr);
		}

    	// set Krylov solver
    	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp_solver_);CHKERRV(ierr);
    	ierr = KSPSetOperators(ksp_solver_,MF_operator_,MF_operator_);CHKERRV(ierr);	    		    	   
    	ierr = KSPSetType(ksp_solver_,ksp_type_);CHKERRV(ierr);     	
    	ierr = KSPSetTolerances(ksp_solver_,cfg.solver.ksp_rtol,PETSC_DEFAULT,PETSC_DEFAULT, cfg.solver.ksp_max_it);CHKERRV(ierr);
    	ierr = KSPSetNormType(ksp_solver_, KSP_NORM_UNPRECONDITIONED);CHKERRV(ierr);    	
    	ierr = KSPGMRESSetRestart(ksp_solver_, cfg.solver.gmres_restart);CHKERRV(ierr);

    	// some warnings 
    	if (mpi_rank_ == 0 and RT_problem_->verbose_)
    	{
    		if (not using_prec_ && std::strcmp(ksp_type_, KSPFGMRES) == 0) std::cout << "WARNING: using FGMRES with no preconditioner, switch to GMRES for better performance." << std::endl;    	
    		if (using_prec_     && std::strcmp(ksp_type_, KSPGMRES)  == 0) std::cout << "WARNING: using GMRES with matrix-free preconditioner can be unsafe, switch to FGMRES." << std::endl;    	
    	}

    	// set preconditioner
    	ierr = KSPGetPC(ksp_solver_,&pc_);CHKERRV(ierr);    		    	    		

    	// set MF_operator_approx_		
    	if (using_prec_)
    	{    	    		
    		if (mf_ctx_.pc_use_J_KQ_)        		
    		{    			
    			ierr = MatCreateShell(PETSC_COMM_WORLD,mf_ctx_.local_J_KQ_size_,mf_ctx_.local_J_KQ_size_,
    												   mf_ctx_.tot_J_KQ_size_,mf_ctx_.tot_J_KQ_size_,
													   (void*)&mf_ctx_,&MF_operator_approx_);CHKERRV(ierr); 
    		}
			else if (mf_ctx_.unpolarized_prec_)
			{				
				ierr = MatCreateShell(PETSC_COMM_WORLD,RT_problem_->local_size_unpolarized_,RT_problem_->local_size_unpolarized_,
													   RT_problem_->tot_size_unpolarized_,RT_problem_->tot_size_unpolarized_,
													   (void*)&mf_ctx_,&MF_operator_approx_);CHKERRV(ierr); 
			}
			else
			{
				ierr = MatCreateShell(PETSC_COMM_WORLD,local_size,local_size,RT_problem_->tot_size_,RT_problem_->tot_size_,
					                                   (void*)&mf_ctx_,&MF_operator_approx_);CHKERRV(ierr); 
			}

			ierr = MatShellSetOperation(MF_operator_approx_,MATOP_MULT,(void(*)(void))UserMult_approx);CHKERRV(ierr);		

			// set PC solver
			ierr = KSPCreate(PETSC_COMM_WORLD,&mf_ctx_.pc_solver_);CHKERRV(ierr);
			// ierr = KSPSetOptionsPrefix(mf_ctx_.pc_solver_, "pc_");CHKERRV(ierr);
    		ierr = KSPSetOperators(mf_ctx_.pc_solver_,MF_operator_approx_,MF_operator_approx_);CHKERRV(ierr);

    		ierr = KSPSetType(mf_ctx_.pc_solver_,pc_ksp_type_);CHKERRV(ierr); 

    		// const PetscInt GMRES_restart = 10;
    		// ierr = KSPGMRESSetRestart(mf_ctx_.pc_solver_, GMRES_restart);CHKERRV(ierr);	

    		// const int max_its = 10;
    		// const Real r_tol = 1e-11;
    		// ierr = KSPSetFromOptions(mf_ctx_.pc_solver_);CHKERRV(ierr);    		
    		ierr = KSPSetTolerances(mf_ctx_.pc_solver_,cfg.prec.pc_rtol,PETSC_DEFAULT,PETSC_DEFAULT, cfg.prec.pc_max_it);CHKERRV(ierr);
    		// ierr = KSPSetTolerances(mf_ctx_.pc_solver_,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT, max_its);CHKERRV(ierr);
    		// ierr = KSPSetNormType(mf_ctx_.pc_solver_, KSP_NORM_NONE);CHKERRV(ierr);
    		// ierr = KSPSetInitialGuessNonzero(mf_ctx_.pc_solver_, PETSC_TRUE);CHKERRV(ierr); // TEST

    		// set PC
    		ierr = PCSetType(pc_,PCSHELL);CHKERRV(ierr);
			ierr = PCShellSetContext(pc_, &mf_ctx_);CHKERRV(ierr);		
			ierr = PCShellSetApply(pc_,MF_pc_Apply);CHKERRV(ierr);				
			ierr = PCShellSetDestroy(pc_,MF_pc_Destroy);CHKERRV(ierr);	

			// adding some options for verbosity
	    	ierr = PetscOptionsSetValue(NULL, "-ksp_converged_reason", "");CHKERRV(ierr);
#define PRECOND_MONITORING // If it dosent work with some petsc versions, switch-off this macro.
#ifdef PRECOND_MONITORING
			ierr = KSPSetFromOptions(mf_ctx_.pc_solver_);
			CHKERRV(ierr);
			if (cfg.prec.verbose)
			{
				ierr = KSPMonitorCancel(mf_ctx_.pc_solver_);
				CHKERRV(ierr);
				ierr = KSPMonitorSet(mf_ctx_.pc_solver_, precond_ksp_monitor, NULL, NULL);
				CHKERRV(ierr);
			} else {
				ierr = KSPMonitorCancel(mf_ctx_.pc_solver_);
				CHKERRV(ierr);
			}			
#endif
		}
    	else
    	{
    		ierr = PCSetType(pc_,PCNONE);CHKERRV(ierr);
    	}


		if (RT_problem_->verbose_)  {
			ierr = PetscOptionsSetValue(NULL, "-ksp_monitor", "");CHKERRV(ierr);
			ierr = PetscOptionsSetValue(NULL, "-ksp_monitor_true_residual", "");CHKERRV(ierr);    // WARNING: this costs
			ierr = PetscOptionsSetValue(NULL, "-ksp_view", "");CHKERRV(ierr);			
		}

		ierr = PetscOptionsSetValue(NULL, "-ksp_converged_reason", "");CHKERRV(ierr);

		// print warning when mponitoring true residual
		PetscBool ksp_true_res_flg;
  		PetscOptionsHasName(NULL,NULL,"-ksp_monitor_true_residual",&ksp_true_res_flg);
        if (ksp_true_res_flg) PetscPrintf(PETSC_COMM_WORLD,"WARNING: -ksp_monitor_true_residual may significantly impact performance.\n\n");

    	// extra options from command line   	
    	ierr = KSPSetFromOptions(ksp_solver_);CHKERRV(ierr);
    	// ierr = KSPSetFromOptions(mf_ctx_.pc_solver_);CHKERRV(ierr);    
    	// ierr = PCSetFromOptions(pc_);CHKERRV(ierr);	     

		// test
		// ierr = MatMult(MF_operator_, rhs_, RT_problem_->I_vec_);CHKERRV(ierr);	        				
		// const std::string filename =  "../output/rhs_" + std::to_string(mpi_size_) + ".m";
  		// const std::string varible  =  "rhs" + std::to_string(mpi_size_);
  		// save_vec(rhs_, filename.c_str(), varible.c_str());        
	}

	// solve linear system
	inline void solve()
	{	
		const double start = MPI_Wtime();	
		
		if (mpi_rank_ == 0) std::cout << "\nStart linear solve..." << std::endl;		

		if (mf_ctx_.ksp_use_J_KQ_)
		{
			solve_with_J_KQ();
		}	
		else
		{			
			PetscErrorCode ierr;			
			ierr = KSPSolve(ksp_solver_, rhs_, RT_problem_->I_vec_);CHKERRV(ierr);
		}

		MPI_Barrier(MPI_COMM_WORLD); 
		if (mpi_rank_ == 0) std::cout << "Solve time (s) = " << MPI_Wtime() - start << std::endl;	

		// update I_field for later use
		mf_ctx_.vec_to_field(RT_problem_->I_field_, RT_problem_->I_vec_);	//TODO	remove
	}

	// To be properly tested 
	inline void solve_with_J_KQ()
	{				
		if (mpi_rank_ == 0) std::cout << "\nStart J_KQ solution: WARNING at the moment this has no eps_csc." << std::endl;			

		PetscErrorCode ierr;
		
		// map rhs to JKQ
		mf_ctx_.I_vec_to_J_KQ_vec(rhs_, mf_ctx_.y_J_KQ_);     

		// solve in JKQ
		ierr = KSPSolve(ksp_solver_, mf_ctx_.y_J_KQ_, mf_ctx_.x_J_KQ_);CHKERRV(ierr); // TOOD

		// (i) compute emissivity
        mf_ctx_.update_emission_J_KQ(mf_ctx_.x_J_KQ_);
        
        // (ii) perform a formal solution        
        if (mpi_rank_ == 0) std::cout << "\nFinal solution to map J_KQ solution to radiation field...";			
        mf_ctx_.formal_solve(RT_problem_->I_field_, RT_problem_->S_field_, 0.0); 

        // map back to y adding rhs for consistency
        mf_ctx_.field_to_vec(RT_problem_->I_field_, RT_problem_->I_vec_);   

        VecAXPY(RT_problem_->I_vec_, 1.0, rhs_);        
	}

	// TEST (test also rhs as initial guess)
	inline void solve_LSA()
	{	
		double start = MPI_Wtime();	

		PetscErrorCode ierr;

		if (mpi_rank_ == 0) std::cout << "\nTesting Last Scattering Approximation..." << std::endl;				
		if (mpi_rank_ == 0) std::cout << "\nStart linear solve..." << std::endl;		
			
		ierr = KSPSolve(ksp_solver_, rhs_, RT_problem_->I_vec_);CHKERRV(ierr);		

		MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank_ == 0) std::cout << "Solve time (s) = " << MPI_Wtime() - start << std::endl;	
		
		// remove polarization
		if (mpi_rank_ == 0) std::cout << "\nRemoving polarization..." << std::endl;	
		set_zero_polarization(RT_problem_->I_vec_);		

		ierr = KSPReset(ksp_solver_);CHKERRV(ierr);    	

		if (mpi_rank_ == 0) std::cout << "Set LSA KSP" << std::endl;
		ierr = KSPCreate(PETSC_COMM_WORLD,&ksp_solver_);CHKERRV(ierr);
    	ierr = KSPSetOperators(ksp_solver_,MF_operator_,MF_operator_);CHKERRV(ierr);	    		    	   
    	ierr = KSPSetType(ksp_solver_,KSPGMRES);CHKERRV(ierr);     	
    	ierr = KSPSetTolerances(ksp_solver_, 1e-6,PETSC_DEFAULT,PETSC_DEFAULT, 10);CHKERRV(ierr);
    	ierr = KSPSetNormType(ksp_solver_, KSP_NORM_UNPRECONDITIONED);CHKERRV(ierr);    	    	
    	ierr = KSPSetInitialGuessNonzero(ksp_solver_, PETSC_TRUE);CHKERRV(ierr); 

		ierr = PetscOptionsSetValue(NULL, "-ksp_monitor", "");CHKERRV(ierr);
        ierr = KSPSetFromOptions(ksp_solver_);CHKERRV(ierr);

		if (mpi_rank_ == 0) std::cout << "Remove prec" << std::endl;	
		ierr = KSPGetPC(ksp_solver_,&pc_);CHKERRV(ierr);    	
		ierr = PCSetType(pc_,PCNONE);CHKERRV(ierr);

		if (mpi_rank_ == 0) std::cout << "Start final Lambdas:" << std::endl;	

		start = MPI_Wtime();	
		ierr = KSPSolve(ksp_solver_, rhs_, RT_problem_->I_vec_);CHKERRV(ierr);

		MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank_ == 0) std::cout << "Solve time LSA (s) = " << MPI_Wtime() - start << std::endl;	

		// update I_field for later use
		mf_ctx_.vec_to_field(RT_problem_->I_field_, RT_problem_->I_vec_);	//TODO	remove?
	}


	// solve linear system with aproxximate emissivity (can be used ofr initial guess)
	inline void solve_approx()
	{	

		// TODO better, with no need for switch... 

		// emissivity for initial guess
		using emission_coefficient_components = rii_include::emission_coefficient_computation::emission_coefficient_components;

	    std::list<emission_coefficient_components> components_approx_IG{    
            emission_coefficient_components::epsilon_R_II_AA_FAST
          , emission_coefficient_components::epsilon_R_III_GL
          , emission_coefficient_components::epsilon_csc
    	};           	
    
    	mf_ctx_.epsilon_fun_approx_ = mf_ctx_.ecc_sh_ptr_->make_computation_function(components_approx_IG);


		const double start = MPI_Wtime();	
		
		if (mpi_rank_ == 0) std::cout << "\nStart approximate linear solve..." << std::endl;		

		PetscErrorCode ierr;			
		ierr = KSPSolve(mf_ctx_.pc_solver_, rhs_, RT_problem_->I_vec_);CHKERRV(ierr);		
		
		MPI_Barrier(MPI_COMM_WORLD); 
		if (mpi_rank_ == 0) std::cout << "Solve time (s) = " << MPI_Wtime() - start << std::endl;			

		// TEST
		// remove polarization
		// if (mpi_rank_ == 0) std::cout << "\nRemoving polarization..." << std::endl;	
		// set_zero_polarization(RT_problem_->I_vec_);		
		ierr = KSPSetInitialGuessNonzero(ksp_solver_, PETSC_TRUE);CHKERRV(ierr); 

		// emissivity for preconidtioner 
	    std::list<emission_coefficient_components> components_approx{emission_coefficient_components::epsilon_pCRD_limit};       
    
    	mf_ctx_.epsilon_fun_approx_ = mf_ctx_.ecc_sh_ptr_->make_computation_function(components_approx);
	}

	
	// set Q,U,V = 0
	inline void set_zero_polarization(Vec& I_vec)
	{
    	PetscErrorCode ierr;
    	PetscScalar* a;
    	PetscInt n_local;

    	ierr = VecGetLocalSize(I_vec, &n_local); CHKERRV(ierr);
    	ierr = VecGetArray(I_vec, &a); CHKERRV(ierr);

    	// Zero Q, U, V in each (I,Q,U,V) block
    	for (PetscInt i = 0; i < n_local; i += 4)
    	{
        	if (i + 1 < n_local) a[i + 1] = 0.0; // Q
        	if (i + 2 < n_local) a[i + 2] = 0.0; // U
        	if (i + 3 < n_local) a[i + 3] = 0.0; // V
    	}

	    ierr = VecRestoreArray(I_vec, &a); CHKERRV(ierr);
	}

	inline void solve_checkpoint(const std::string output_path, const int checkpoint_interval)
	{
		PetscErrorCode ierr;
		KSPConvergedReason reason = KSP_DIVERGED_ITS;
		std::string output_file;
    		
		int counter = 0;
		PetscInt its;

		ierr = KSPSetInitialGuessNonzero(ksp_solver_, PETSC_TRUE);CHKERRV(ierr); 
		ierr = KSPSetTolerances(ksp_solver_,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT, checkpoint_interval);CHKERRV(ierr);						

		// while KSPSolve reaches max it
		while (reason == KSP_DIVERGED_ITS)
		{
			ierr = KSPSolve(ksp_solver_, rhs_, RT_problem_->I_vec_);CHKERRV(ierr); // TODO RT_problem_->I_vec_ => field.getVec()
			ierr = KSPGetConvergedReason(ksp_solver_, &reason);CHKERRV(ierr); 
			ierr = KSPGetIterationNumber(ksp_solver_, &its);CHKERRV(ierr); 

			counter += its;

			// update I_field for write_surface_point_profiles()
			mf_ctx_.vec_to_field(RT_problem_->I_field_, RT_problem_->I_vec_); // TODO remove	
			
			output_file = output_path + "CP" + std::to_string(counter);
			RT_problem_->write_surface_point_profiles(output_file, 0, 0);

		}		
		
		// update I_field for later use
		mf_ctx_.vec_to_field(RT_problem_->I_field_, RT_problem_->I_vec_);	//TODO	remove
	}

	//TODO 
	inline void apply_formal_solver()
	{		
		double start = MPI_Wtime();		

		// // set source fun
		// if (mpi_rank_ == 0) std::cout << "WARNING: setting source function in apply_formal_solver()" << std::endl;

		// auto S_dev = RT_problem_->S_field_;

		// RT_problem_->space_grid_->parallel_for([&](int i, int j, int k) 
		// {         
		// 	auto block = S_dev->block(i, j, k);

		// 	for (int b = 0; b < (int)RT_problem_->block_size_; ++b) 
		// 	{	    
		// 		block[b] = 0.1;     	   	
		// 		// if (block[b] != 0) std::cout << "S not zero!" << std::endl;
		// 	}
		// });

		if (mpi_rank_ == 0) std::cout << "Start formal solve..." << std::endl;
		
		mf_ctx_.formal_solve_global(RT_problem_->I_field_, RT_problem_->S_field_, 1.0);		
						
		MPI_Barrier(MPI_COMM_WORLD); double end = MPI_Wtime();
		if (mpi_rank_ == 0) std::cout << "Formal solve time (s) = " << end - start << std::endl;	

		// update I_vec for later use
		mf_ctx_.field_to_vec(RT_problem_->I_field_, RT_problem_->I_vec_); // TODO remove
	}	

	//TODO
	inline void compute_emission()
	{
		double start = MPI_Wtime();		
		if (mpi_rank_ == 0) std::cout << "Computing emission..." << std::endl;

		// // test
  		// VecSet(RT_problem_->I_vec_,0.0);
		
		// compute new emission in S_field_ 
  		mf_ctx_.update_emission(RT_problem_->I_vec_);   
    		
		MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank_ == 0) std::cout << "Computing emission took (s) = " << MPI_Wtime() - start << std::endl;	

		// start = MPI_Wtime();		

		// mf_ctx_.update_emission(RT_problem_->I_vec_, true); 

		// MPI_Barrier(MPI_COMM_WORLD);
		// if (mpi_rank_ == 0) std::cout << "Computing approximate emission took (s) = " << MPI_Wtime() - start << std::endl;	

		// // use I_vec to store S for later use
		// mf_ctx_.field_to_vec(RT_problem_->S_field_, RT_problem_->I_vec_);
	}


	//TODO
	// set the radiation field in an arbitrary direction Omega
	inline void apply_formal_solver_Omega(const double theta, const double chi)
	{
		// allocate new data structure
		if (not mf_ctx_.formal_solution_Omega_)
		{
			if (mpi_rank_ == 0) std::cout << "\nAllocating fields for new direction...";

			RT_problem_->allocate_fields_Omega();
			mf_ctx_.init_serial_fields_Omega();
			mf_ctx_.formal_solution_Omega_ = true;

			if (mpi_rank_ == 0) std::cout << "done" << std::endl;
		}		

		// set eta and rhos 
	    RT_problem_->set_eta_and_rhos_Omega(theta, chi);
	
		const double clock_start = MPI_Wtime();				

		// update emissivity with current I_field (in all directions)
		mf_ctx_.update_emission_Omega(RT_problem_->I_vec_, theta, chi);	// TODO RT_problem_->I_vec_ => field->getVec()	

		const double clock_end = MPI_Wtime();
		const double clock_diff = clock_end - clock_start;

		if (mpi_rank_ == 0){ 
			std::cout << "\nComputing emission took (s) = " << clock_diff << "    file: " << __FILE__ << ":" << __LINE__ << std::endl;
			std::cout.flush();
		}
		
		// formal solve
		MPI_Barrier(MPI_COMM_WORLD);

		if (mpi_rank_ == 0) std::cout << "Start formal solve in Omega..." << "    file: " << __FILE__ << ":" << __LINE__ << std::endl;
		
		double clock_start_formal = MPI_Wtime();
		mf_ctx_.formal_solve_ray(theta, chi);

		MPI_Barrier(MPI_COMM_WORLD);
		double clock_end_formal = MPI_Wtime();

		double clock_diff_formal = clock_end_formal - clock_start_formal;
		if (mpi_rank_ == 0) std::cout << "Formal solve time (s) = " << clock_diff_formal << "    file: " << __FILE__ << ":" << __LINE__ << std::endl;
	}


	inline void get_formation_height(const double theta, const double chi)
	{
		mf_ctx_.get_formation_height(theta, chi);
	}

	// inline void test()
	// {		
	// 	const auto N_nu  = RT_problem_->N_nu_;
	// 	const auto N_pol = RT_problem_->N_pol_;
	// 	const auto space_grid = RT_problem_->space_grid_;
	// 	const auto space_grid_s = mf_ctx_.space_grid_serial_;

	// 	const int block_size = N_pol * N_nu;
	// 	const int local_block_size = N_pol;

    //     Field_ptr_t test_field_        = std::make_shared<Field>("test"       , space_grid,   N_pol, std::vector<PetscInt>{N_nu}, false);       
    //     Field_ptr_t test_field_serial_ = std::make_shared<Field>("test_serial", space_grid_s, local_block_size, false);
               
    //     ReMap3D test_remap;
    //     test_remap.init(space_grid, space_grid_s, block_size, local_block_size);         
       
    // 	// write eta to the serial grid
    // 	if (mpi_rank_ == 0) std::cout << "Sending eta to serial" << std::endl;
    // 	MPI_Barrier(MPI_COMM_WORLD);
    // 	test_remap.from_space_to_block_distributed(test_field_, test_field_serial_);         
    // 	MPI_Barrier(MPI_COMM_WORLD);
    // 	if (mpi_rank_ == 0) std::cout << "DONE" << std::endl;
	// }
	

	inline void test_transfer()
	{			
		PetscErrorCode ierr; 

		const auto g_dev      = RT_problem_->space_grid_;
		const auto field_dev  = RT_problem_->I_field_;	
		const auto block_size = RT_problem_->block_size_;

		// indeces
		const auto [i_start, j_start, k_start] = g_dev->getGhostMargins(); 

		const int i_end = i_start + g_dev->getLocalSizeX();
		const int j_end = j_start + g_dev->getLocalSizeY();
		const int k_end = k_start + g_dev->getLocalSizeZ();

		int counter = 0;

		PetscInt istart;	
		ierr = VecGetOwnershipRange(RT_problem_->I_vec_, &istart, NULL);CHKERRV(ierr);	

		// init
		for (int k = k_start; k < k_end; ++k)					
		{															
			for (int j = j_start; j < j_end; ++j)
			{
				for (int i = i_start; i < i_end; ++i)				
				{
					for (int b = 0; b < block_size; b++) 
					{			
						field_dev->block(i, j, k)[b] = istart + counter;							

						counter++;
					}							
				}
			}
		}

		// test 1
		mf_ctx_.field_to_vec(RT_problem_->I_field_, RT_problem_->I_vec_); // TODO remove this

		save_vec(RT_problem_->I_vec_,  "../output/vec1.m" ,"vec_1");   

		// test 2 
		mf_ctx_.vec_to_field(RT_problem_->I_field_, RT_problem_->I_vec_); // TODO remove this 
		mf_ctx_.field_to_vec(RT_problem_->I_field_, RT_problem_->I_vec_); // TODO remove this

		save_vec(RT_problem_->I_vec_,  "../output/vec2.m" ,"vec_2");   

		// const std::string filename =  "../output/I_" + std::to_string(mpi_size_) + ".m";
	}

	inline void free_fields_memory()
	{
		if (mpi_rank_ == 0) std::cout << "Freeing RT_solver fields memory..." << std::endl;				

		mf_ctx_.I_field_serial_.reset();
		mf_ctx_.S_field_serial_.reset();
		mf_ctx_.eta_field_serial_.reset();
		mf_ctx_.rho_field_serial_.reset();

		PetscErrorCode ierr;

		if (mf_ctx_.unpolarized_prec_)
		{			
			mf_ctx_.I_unpol_field_serial_.reset();
			mf_ctx_.S_unpol_field_serial_.reset();		

			RT_problem_->I_unpol_field_.reset();
    		RT_problem_->S_unpol_field_.reset();

			ierr = VecDestroy(&(mf_ctx_.x_unpol_));CHKERRV(ierr);
			ierr = VecDestroy(&(mf_ctx_.y_unpol_));CHKERRV(ierr);
			ierr = VecDestroy(&(mf_ctx_.x_pol_));CHKERRV(ierr);
			ierr = VecDestroy(&rhs_);CHKERRV(ierr);
		}
		
		// ierr = MatDestroy(&MF_operator_);CHKERRV(ierr);
		// ierr = MatDestroy(&MF_operator_approx_);CHKERRV(ierr);
		// ierr = KSPDestroy(&ksp_solver_);CHKERRV(ierr);
		// ierr = PCDestroy(&pc_);CHKERRV(ierr);
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
	KSPType ksp_type_;    //= KSPFGMRES; //KSPFBCGS //KSPRICHARDSON; //KSPPIPEFGMRES
	KSPType pc_ksp_type_; //= KSPGMRES;  //KSPPIPEFGMRES //KSPBCGS
	PC pc_;
	
	bool using_prec_;	
			
	// assemble Lam[eps_th] + t
	void assemble_rhs();		

	void print_info();	
};


#endif 
// 
