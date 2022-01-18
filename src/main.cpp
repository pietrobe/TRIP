#include "RT_solver.hpp"
#include <chrono>  

//  make -j4 && mpirun -n 4 ./main

int main(int argc, char *argv[])
{    	
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    {	    
	    const size_t N_theta = 2; 
	    const size_t N_chi   = 2;     
	    
	    RT_problem rt_problem("../input/FAL-C/2_B0_V10P_12T_8C_99F", N_theta, N_chi);	

	    auto rt_problem_ptr = std::make_shared<RT_problem>(rt_problem);

	    RT_solver rt_solver(rt_problem_ptr, "DELO_linear");
	    rt_solver.solve();
	}
    
	Kokkos::finalize();
    return MPI_Finalize();
}


// output 
// python ../../sgrid/scripts/transpose_data.py -x 4 -y 4 -z 70 -b 99 -p sigma.raw 