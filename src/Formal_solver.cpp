#include "Formal_solver.hpp"

// solve I' = KI - S

void Formal_solver::one_step(const Real dt, input_vec &K1, input_vec &K2, input_vec &S1, input_vec &S2, input_vec &I_in, output_vec &I_out){

	if (debug_mode_ and (K2.size() != 16 or S2.size() != 4 or I_in.size() != 4) and mpi_rank_ == 0) std::cerr << "ERROR: wrong input size in one_step().\n";
	
	std::vector<Real> A(16); 
	std::vector<Real> b(4); 
	
	if (type_ == "implicit_Euler") // K1 and S1 are not needed here 
	{		
		// compute b = I_in - dt * S and A = Id - dt * K		
		int index_ij;

		for (int i = 0; i < 4; ++i) 
		{
			b[i] = I_in[i] - dt * S2[i];			

			for (int j = 0; j < 4; ++j)
			{
				index_ij = 4 * i + j;

				A[index_ij] = (i == j) - dt * K2[index_ij];
			}			
		}
			
		I_out = solve_4_by_4_system(A, b);			
	}
	else if (type_ == "trapezoidal" or type_ == "Crank–Nicolson")
	{		
		if (debug_mode_ and (K1.size() != 16 or S1.size() != 4) and mpi_rank_ == 0) std::cerr << "\nERROR: wrong input size for Crank–Nicolson!\n";

		std::vector<Real> K_times_I_in(4); 

		Real dt_half = 0.5 * dt;

		int index_ij;

		for (int i = 0; i < 4; ++i) 			
		{						
			for (int j = 0; j < 4; ++j)
			{
				index_ij = 4 * i + j;

				// Id - (dt/2) * K
				A[index_ij] = (i == j) - dt_half * K2[index_ij];

				// K * I_in
				K_times_I_in[i] += K1[index_ij] * I_in[j]; 
			}	

			b[i] = I_in[i] + dt_half * (K_times_I_in[i] - S1[i] - S2[i]);			
		}
			
		I_out = solve_4_by_4_system(A, b);	
	}
	else if (type_ == "DELO_linear")
	{    			
		// solve I' = - KI + S
		const long double dt_aux = -dt;         
        const long double E = std::exp(-dt_aux);      
        const long double F = 1.0 - E;
        const long double G = (1.0 - (1.0 + dt_aux) * E) / dt_aux; 
        
        std::vector<Real> K_times_I_in(4, 0.0);
        int index_ij;
				
        bool Id;

     	for (int i = 0; i < 4; ++i) 
		{			
			for (int j = 0; j < 4; ++j)
			{
				index_ij = 4 * i + j;

				Id = (i == j);
			
				A[index_ij] = Id + (F - G) * (K2[index_ij] - Id);

				// K * I_in
				K_times_I_in[i] += (K1[index_ij] - Id) * I_in[j]; 								
			}			

			b[i] = E * I_in[i] - G * K_times_I_in[i] + G * S1[i] + (F - G) * S2[i];
		}

		I_out = solve_4_by_4_system(A, b);	
		// I_out = solve_4_by_4_system_optimized(A,b);
	}
	else	
	{
		if (mpi_rank_ == 0) std::cerr << "\nERROR: " << type_ << " is not supported as formal solver!\n";
	}

	// check valore negativo
	if (debug_mode_ and I_out[0] < 0 and I_in[0] >= 0 and S1[0] > 0  and S2[0] > 0) std::cout << "\nWARNING: negative intensity detected!\n";	

	// if (I_out[0] < 0) I_out[0] = 0;		
}



void Formal_solver::solve(input_vec &dts, input_field &K, input_field &S, input_vec &I_in, output_field &I_out){

	size_t N_tau = dts.size() + 1;

	// check input sizes
	if (debug_mode_ and (K.size() != N_tau or S.size() != N_tau or I_in.size() != 4) and mpi_rank_ == 0) std::cerr << "\nERROR: wrong input size in Formal_solver::solve()!\n";

	I_out.resize(N_tau);

	// set initial condition 
	I_out[0] = I_in;

	// solve ODE
	for (size_t i = 1; i < N_tau; ++i)
	{
		one_step(dts[i-1], K[i-1], K[i], S[i-1], S[i], I_out[i-1], I_out[i]);		
	}
}











