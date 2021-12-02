#ifndef Formal_solver_hpp
#define Formal_solver_hpp

#include "Utilities.hpp"

typedef const std::vector<double> input_vec;
typedef const std::vector<std::vector<double> > input_field;

typedef std::vector<double> output_vec;
typedef std::vector<std::vector<double> > output_field;

// base class
class Formal_solver
{
public:

	Formal_solver() = default;

	Formal_solver(std::string type, bool debug_mode = true) : type_(type), debug_mode_(debug_mode)
	{				
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);		
		
		// print type and check
		if (mpi_rank_ == 0)
		{
			if (type_ == "implicit_Euler")
			{
				std::cout << "Formal solver: implicit Euler\n";	
			}
			else if (type_ == "trapezoidal" or type_ == "Crank–Nicolson")
			{
				std::cout << "Formal solver: Crank–Nicolson\n";	
			}
			else if (type_ == "DELO_linear")
			{
				std::cout << "Formal solver: DELO linear\n";	
			}
			else
			{
				std::cerr << "ERROR: " << type_ << " is not supported as formal solver.\n";
			}

			if (debug_mode_) std::cout << "Formal solver bebug mode enabled.\n";	
		}				
	}

	// solve on a ray
	void solve(input_vec &dts, input_field &K, input_field &S, input_vec &I_in, output_field &I_out);

	// solve, for one step (dt), I' = K * I - S with initial condition I_in and K = [K1 K2] and S = [S1 S2]
	void one_step(const double dt, input_vec &K1, input_vec &K2, input_vec &S1, input_vec &S2, input_vec &I_in, output_vec &I_out);

private:

	// MPI variables
	int mpi_rank_;

	// type of formal solver
	std::string type_;

	// debug flag
	bool debug_mode_;
};
	
#endif 