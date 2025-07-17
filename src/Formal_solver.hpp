#ifndef Formal_solver_hpp
#define Formal_solver_hpp

#include "Utilities.hpp"

typedef const std::vector<Real> input_vec;
typedef const std::vector<std::vector<Real> > input_field;

typedef std::vector<Real> output_vec;
typedef std::vector<std::vector<Real> > output_field;


// base class
class Formal_solver
{
public:

	Formal_solver() = default;

	Formal_solver(std::string type, bool debug_mode = true) : type_(type), debug_mode_(debug_mode)
	{				
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);		
		
		// print type and check		
		stencil_size_ = 2;

		if (type_ == "implicit_Euler")
		{
			if (mpi_rank_ == 0) std::cout << "Formal solver: implicit Euler\n";	
		}
		else if (type_ == "trapezoidal" or type_ == "Crank–Nicolson")
		{
			if (mpi_rank_ == 0) std::cout << "Formal solver: Crank–Nicolson\n";	
		}
		else if (type_ == "DELO_linear")
		{
			if (mpi_rank_ == 0) std::cout << "Formal solver: DELO linear\n";	
		}
		else if (type_ == "SC_linear")
		{
			if (mpi_rank_ == 0) std::cout << "Unpolarized formal solver: linear SC\n";	
		}
		else if (type_ == "SC_parabolic")
		{
			if (mpi_rank_ == 0) std::cout << "Unpolarized formal solver: parabolic SC\n";	

			stencil_size_ = 3;
		}
		else if (type_ == "BESSER")
		{
			if (mpi_rank_ == 0) std::cout << "Formal solver: BESSER\n";	

			stencil_size_ = 3;
		}
		else
		{
			if (mpi_rank_ == 0) std::cerr << "ERROR: " << type_ << " is not supported as formal solver.\n";
			if (mpi_rank_ == 0) std::cerr << "Supported inputs: implicit_Euler, Crank–Nicolson, DELO_linear, SC_linear, SC_parabolic, and BESSER.\n";
		}

		if (debug_mode_ and mpi_rank_ == 0) std::cout << "Formal solver bebug mode enabled.\n";				
	}

	// getters
	inline int         get_stencil_size() const {return stencil_size_;}
	inline std::string get_type()         const {return type_;        }

	// solve on a ray
	void solve(input_vec &dts, input_field &K, input_field &S, input_vec &I_in, output_field &I_out);

	// solve, for one step (dt), I' = K * I - S with initial condition I_in and K = [K1 K2] and S = [S1 S2]
	void one_step(const Real dt, input_vec &K1, input_vec &K2, input_vec &S1, input_vec &S2, input_vec &I_in, output_vec &I_out);


	// one step method for quadratic stencils
	void one_step_quadratic(const Real dt_1, const Real dt_2, input_vec &K1, input_vec &K2, input_vec &K3,
								 						      input_vec &S1, input_vec &S2, input_vec &S3,
								 							  input_vec &I_in, output_vec &I_out);

	// in the case of scalar equations (just intensity)
	Real one_step(const Real dt, const Real I_in, const Real S1, const Real S2);
	Real one_step_quadratic(const Real dt1, const Real dt2, const Real I_in, const Real S1, const Real S2, const Real S3);

private:

	// MPI variables
	int mpi_rank_;

	// stencil size: 2 for linear, 3 for quadratic
	int stencil_size_;

	// type of formal solver
	std::string type_;

	// debug flag
	bool debug_mode_;
};
	
#endif 