#include "RT_problem.hpp"
#include <chrono>  

//  make -j4 && mpirun -n 4 ./main

int main(int argc, char *argv[])
{    	
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    {
	    int mpi_rank;

	    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	    const size_t N_theta = 4; 
	    const size_t N_chi   = 4;     
	    
	    RT_problem rt_problem("../input/FAL-C/1_B0_V0_12T_8C_99F", N_theta, N_chi);	    
	}
    
	Kokkos::finalize();
    return MPI_Finalize();
}
