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

	    const size_t N_theta = 10; 
	    const size_t N_chi   = 10;     

	    const std::string input_path= "path/to/input/files";

	    RT_problem rt_problem(input_path, N_theta, N_chi);

	    if (mpi_rank == 0) std::cout << "RT problem constructed" << std::endl;		
	}
    
	Kokkos::finalize();
    return MPI_Finalize();
}
