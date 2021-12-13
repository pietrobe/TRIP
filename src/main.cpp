#include "RT_problem.hpp"
#include <chrono>  

//  make -j4 && mpirun -n 4 ./main

int main(int argc, char *argv[])
{    	
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    {	    
	    const size_t N_theta = 2; 
	    const size_t N_chi   = 2;     
	    
	    RT_problem rt_problem("../input/FAL-C/1_B0_V0_12T_8C_99F", N_theta, N_chi);				   	
	}
    
	Kokkos::finalize();
    return MPI_Finalize();
}


// output 
// python ../../sgrid/scripts/transpose_data.py -x 4 -y 4 -z 70 -b 99 -p sigma.raw 