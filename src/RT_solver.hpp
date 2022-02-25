#ifndef RT_solver_hpp
#define RT_solver_hpp

#include "Formal_solver.hpp"
#include "RT_problem.hpp"
// #include <rii_emission_coefficient.h>
// #include <thread>

extern PetscErrorCode UserMult(Mat mat,Vec x,Vec y);
extern PetscErrorCode UserMult_approx(Mat mat,Vec x,Vec y);
// extern PetscErrorCode MF_pc_Destroy(PC pc);
// extern PetscErrorCode MF_pc_Apply(PC pc,Vec x,Vec y);

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
	
	// // pointer for emission module and offset
	// std::shared_ptr<rii_include::emission_coefficient_computation> ecc_sh_ptr_;
	// rii_include::offset_function_cartesian offset_f_;	
	// rii_include::emission_coefficient_computation::compute_height_function_type epsilon_computation_function_;
	// rii_include::emission_coefficient_computation::compute_height_function_type epsilon_computation_function_approx_;

	// change data format
	void field_to_vec(const Field_ptr_t field, Vec &v);
	void vec_to_field(Field_ptr_t field, const Vec &v);
	
	// find intersection
	void find_intersection(double theta, double chi, const double Z_down, const double Z_top, const double L, t_intersect *T);
	std::vector<t_intersect> find_prolongation(double theta, double chi, const double dz, const double L);
	std::vector<double> long_ray_steps(const std::vector<t_intersect> T, const Field_ptr_t I_field, const Field_ptr_t S_field, const int i, const int j, const int k, const int block_index);

	void get_2D_weigths(const double x, const double y, double *w);

	// formal solvers methods 	
	void formal_solve_local( Field_ptr_t I_field, const Field_ptr_t S_field, const Real I0);
	void formal_solve_global(Field_ptr_t I_field, const Field_ptr_t S_field, const Real I0);		
	
	void apply_bc(Field_ptr_t I_field, const Real I0);	
		
	// emission module from Simone
	void set_up_emission_module();
		
	// update emission in all spatial points (given the current I_field_, update S_field_)
	void update_emission(const Vec &I_field, const bool approx = false);
};

class RT_solver
{
public:
	RT_solver(const std::shared_ptr<RT_problem> RT_problem, input_string formal_solver = "implicit_Euler", const bool using_prec = true, const bool threaded = false) 
	{		
		// assign MPI varaibles and init mf_ctx_
    	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
    	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);  

    	RT_problem_ = RT_problem;  
    	using_prec_ = using_prec;    	

    	mf_ctx_.RT_problem_    = RT_problem;  
    	mf_ctx_.mpi_rank_      = mpi_rank_;
    	mf_ctx_.mpi_size_      = mpi_size_;    	
    	mf_ctx_.formal_solver_ = Formal_solver(formal_solver);     
    	mf_ctx_.set_up_emission_module();  	    		

    	// assemble rhs
    	assemble_rhs();
    	// VecView(rhs_, PETSC_VIEWER_STDOUT_SELF);
  
    	// set linear system
    	PetscErrorCode ierr;

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
			// // set MF_operator_approx_		
			// ierr = MatCreateShell(PETSC_COMM_WORLD,local_size,local_size,RT_problem_->tot_size_,RT_problem_->tot_size_,(void*)&mf_ctx_,&MF_operator_approx_);CHKERRV(ierr); 
			// ierr = MatShellSetOperation(MF_operator_approx_,MATOP_MULT,(void(*)(void))UserMult_approx);CHKERRV(ierr);		

			// // set PC solver 
			// ierr = KSPCreate(PETSC_COMM_WORLD,&mf_ctx_.pc_solver_);CHKERRV(ierr);
   //  		ierr = KSPSetOperators(mf_ctx_.pc_solver_,MF_operator_approx_,MF_operator_approx_);CHKERRV(ierr);	    		
   //  		ierr = KSPSetType(mf_ctx_.pc_solver_,KSPGMRES);CHKERRV(ierr); 

   //  		// ierr = KSPSetFromOptions(mf_ctx_.pc_solver_);CHKERRV(ierr);
   //  		// ierr = KSPSetTolerances(mf_ctx_.pc_solver_,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRV(ierr);

   //  		// set PC
   //  		ierr = PCSetType(pc_,PCSHELL);CHKERRV(ierr);
			// ierr = PCShellSetContext(pc_, &mf_ctx_);CHKERRV(ierr);		
			// ierr = PCShellSetApply(pc_,MF_pc_Apply);CHKERRV(ierr);				
			// ierr = PCShellSetDestroy(pc_,MF_pc_Destroy);CHKERRV(ierr);	
    	}
    	else
    	{
    		ierr = PCSetType(pc_,PCNONE);CHKERRV(ierr);
    	}
    	
    	// extra options from command line   	
    	ierr = KSPSetFromOptions(ksp_solver_);CHKERRV(ierr);
    	ierr = PCSetFromOptions(pc_);CHKERRV(ierr);	    		    
	}

	// solve linear system
	inline void solve()
	{	
		PetscErrorCode ierr;

		Real start, end;
		
		start = MPI_Wtime();							

		// if (mpi_rank_ == 0) cout << "\nStart linear solve..." << endl;	
		// // ierr = KSPSetInitialGuessNonzero(ksp_solver_, PETSC_TRUE);CHKERRV(ierr);	/// -------> TODO check
		// ierr = KSPSolve(ksp_solver_, rhs_, RT_problem_->I_vec_);CHKERRV(ierr);

		// MPI_Barrier(MPI_COMM_WORLD); end = MPI_Wtime();
		// if (mpi_rank_ == 0) cout << "Solve time (s) = " << end - start << endl;

		/////////////////////////////////////////////////////////////
		// test formal solver
		mf_ctx_.apply_bc(RT_problem_->I_field_, 1.0);	

		// RT_problem_->I_field_->write("I_in.raw");		

		// ierr = VecSet(RT_problem_->S_vec_, 1.0);CHKERRV(ierr);
		// mf_ctx_.vec_to_field(RT_problem_->S_field_, RT_problem_->S_vec_); 		
				
		if (mpi_rank_ == 0) cout << "Global formal solve " << endl;
		mf_ctx_.formal_solve_global(RT_problem_->I_field_, RT_problem_->S_field_, 1.0);		
			
		MPI_Barrier(MPI_COMM_WORLD); end = MPI_Wtime();
		if (mpi_rank_ == 0) cout << "Solve time (s) = " << end - start << endl;
	
		// // RT_problem_->I_field_->write("I_out.raw");			
	}
	
private:	

	// MPI varables
	int mpi_rank_;
	int mpi_size_;

	std::shared_ptr<RT_problem> RT_problem_;	

	// MF context
	MF_context mf_ctx_;
	
	// linear system quantities
	Mat MF_operator_;
	// Mat MF_operator_approx_;
	Vec rhs_;
	
	KSP ksp_solver_;
	KSPType ksp_type_ = KSPGMRES;
	PC pc_;

	bool LTE_        = false;
	bool using_prec_ = false;
			
	// assemble Lam[eps_th] + t
	void assemble_rhs();		
};


#endif 
// 