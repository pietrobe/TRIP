#ifndef RT_solver_hpp
#define RT_solver_hpp

#include "Formal_solver.hpp"
#include "RT_problem.hpp"
// #include <rii_emission_coefficient.h>
// #include <thread>

// extern PetscErrorCode UserMult(Mat mat,Vec x,Vec y);
// extern PetscErrorCode UserMult_approx(Mat mat,Vec x,Vec y);
// extern PetscErrorCode MF_pc_Destroy(PC pc);
// extern PetscErrorCode MF_pc_Apply(PC pc,Vec x,Vec y);

// struct for ray - grid intersection 
typedef struct t_intersect {
    int ix[4], iy[4], iz[4];
    double w[4];
    double distance;
} t_intersect;


// matrix-free (MF) structure
struct MF_context {

	std::shared_ptr<RT_problem> RT_problem_;	

	Formal_solver formal_solver_;

	// preconditioner data structures 
	KSP pc_solver_;

	// MPI varables
	int mpi_rank_;
	int mpi_size_;

	bool threaded_;

	// // pointer for emission module and offset
	// std::shared_ptr<rii_include::emission_coefficient_computation> ecc_sh_ptr_;
	// rii_include::offset_function_cartesian offset_f_;	
	// rii_include::emission_coefficient_computation::compute_height_function_type epsilon_computation_function_;
	// rii_include::emission_coefficient_computation::compute_height_function_type epsilon_computation_function_approx_;

	// find intersection
	void find_intersection(double theta, double chi, const double Z_down, const double Z_top, const double L, t_intersect *T);

	// formal solvers methods 	
	void formal_solve_local(Field_ptr_t I_field, const Field_ptr_t S_field, const Real I0);	
	// void formal_solve(Vec &I_field, const Vec &S_field, const double I0);	
	
	void apply_bc(Field_ptr_t I_field, const Real I0);	

	// // TODO
	// void formal_solve_3D(Vec &I_field, const Vec &S_field, const double I0); 

	// // if I < 0 set I = 0	
	// void make_intensity_positive(Vec &I_field);
	
	// // emission module from Simone
	// void set_up_emission_module();
		
	// // update emission in all spatial points (given the current I_field_, update S_field_)
	// void update_emission(const Vec &I_field, const bool approx = false);
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
    	mf_ctx_.threaded_      = threaded;
    	mf_ctx_.formal_solver_ = Formal_solver(formal_solver);     
    	// mf_ctx_.set_up_emission_module();  	    		

  //   	// assemble rhs
  //   	assemble_rhs();
  //   	// save_vec(rhs_, "../output/rhs.m" ,"b");      

  		// mf_ctx_.formal_solve_local(RT_problem_->I_field_, RT_problem_->S_field_, 1.0);	
  		// RT_problem_->I_field_->write("I.raw");		

  //   	// set linear system
  //   	PetscErrorCode ierr;

		// int local_size;
		// ierr = VecGetLocalSize(rhs_, &local_size);CHKERRV(ierr); 		

		// // init user defined Mat mult
		// ierr = MatCreateShell(PETSC_COMM_WORLD,local_size,local_size,RT_problem_->tot_size_,RT_problem_->tot_size_,(void*)&mf_ctx_,&MF_operator_);CHKERRV(ierr); 
		// ierr = MatShellSetOperation(MF_operator_,MATOP_MULT,(void(*)(void))UserMult);CHKERRV(ierr);

  //   	// set Krylov solver
  //   	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp_solver_);CHKERRV(ierr);
  //   	ierr = KSPSetOperators(ksp_solver_,MF_operator_,MF_operator_);CHKERRV(ierr);	    		
  //   	ierr = KSPSetType(ksp_solver_,KSPGMRES);CHKERRV(ierr); 
  //   	// ierr = KSPSetType(ksp_solver_,KSPRICHARDSON);CHKERRV(ierr); // Lambda iteration

  //   	// set preconditioner
  //   	ierr = KSPGetPC(ksp_solver_,&pc_);CHKERRV(ierr);    		    	    		

  //   	if (using_prec_)
  //   	{    	    		
		// 	// set MF_operator_approx_		
		// 	ierr = MatCreateShell(PETSC_COMM_WORLD,local_size,local_size,RT_problem_->tot_size_,RT_problem_->tot_size_,(void*)&mf_ctx_,&MF_operator_approx_);CHKERRV(ierr); 
		// 	ierr = MatShellSetOperation(MF_operator_approx_,MATOP_MULT,(void(*)(void))UserMult_approx);CHKERRV(ierr);		

		// 	// set PC solver 
		// 	ierr = KSPCreate(PETSC_COMM_WORLD,&mf_ctx_.pc_solver_);CHKERRV(ierr);
  //   		ierr = KSPSetOperators(mf_ctx_.pc_solver_,MF_operator_approx_,MF_operator_approx_);CHKERRV(ierr);	    		
  //   		ierr = KSPSetType(mf_ctx_.pc_solver_,KSPGMRES);CHKERRV(ierr); 

  //   		// ierr = KSPSetFromOptions(mf_ctx_.pc_solver_);CHKERRV(ierr);
  //   		// ierr = KSPSetTolerances(mf_ctx_.pc_solver_,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRV(ierr);

  //   		// set PC
  //   		ierr = PCSetType(pc_,PCSHELL);CHKERRV(ierr);
		// 	ierr = PCShellSetContext(pc_, &mf_ctx_);CHKERRV(ierr);		
		// 	ierr = PCShellSetApply(pc_,MF_pc_Apply);CHKERRV(ierr);				
		// 	ierr = PCShellSetDestroy(pc_,MF_pc_Destroy);CHKERRV(ierr);	
  //   	}
  //   	else
  //   	{
  //   		ierr = PCSetType(pc_,PCNONE);CHKERRV(ierr);
  //   	}
    	
  //   	// extra options    	
  //   	ierr = PCSetFromOptions(pc_);CHKERRV(ierr);	
  //   	ierr = KSPSetFromOptions(ksp_solver_);CHKERRV(ierr);	    	
	}

	// solve linear system
	inline void solve()
	{				
		mf_ctx_.apply_bc(RT_problem_->I_field_, 1.0);	

		RT_problem_->I_field_->write("I_in.raw");		

		const int n_iter = 1;		
		
		for (int i = 0; i < n_iter; ++i)
		{
			cout << "local formal solve " << i << endl;
			mf_ctx_.formal_solve_local(RT_problem_->I_field_, RT_problem_->S_field_, 1.0);	
		}	

		RT_problem_->I_field_->write("I_out.raw");		
	}

	// void set_I_from_input(const std::string input_path, Vec &I);

	// // save matrices for Matlab
	// void save_Lamda();
	
private:	

	// MPI varables
	int mpi_rank_;
	int mpi_size_;

	std::shared_ptr<RT_problem> RT_problem_;	
	
	// MF context
	MF_context mf_ctx_;

	// linear system quantities
	Mat MF_operator_;
	Mat MF_operator_approx_;
	Vec rhs_;
	
	KSP ksp_solver_;
	PC pc_;

	bool LTE_ = false;
	bool using_prec_;
			
	// assemble Lam[eps_th] + t
	void assemble_rhs();		
};


#endif 
// 