#include "RT_solver.hpp"
#include "Test_rii_include.hpp"
#include <chrono>
#include "tools.h"

// compile and run with:
// make -j 32 && srun -n 4 ./main 2 4

int main(int argc, char *argv[]) {

  MPI_CHECK(MPI_Init(&argc, &argv));
  PetscInitialize(&argc, &argv, (char *)0, NULL);
  Kokkos::initialize(argc, argv);

  {
    const bool save_output = false;

    if (argc < 2) std::cout << "Missing arguments, usage: srun -n 1 ./solar_3D N_theta N_chi" << std::endl;

    const size_t N_theta = atoi(argv[1]);
    const size_t N_chi   = atoi(argv[2]);

    // const std::string input_path  = "../input/FAL-C/1_B0_V0_12T_8C_99F";
    const std::string input_path  = "../input/FAL-C/1_B0_V0_12T_8C_64F_64";
    // const std::string input_path  = "../input/FAL-C/1_B0_V0_12T_8C_4F_64";

    auto rt_problem_ptr = std::make_shared<RT_problem>(input_path, N_theta, N_chi);

    RT_solver rt_solver(rt_problem_ptr, "DELO_linear");    

    rt_solver.solve();  

    rt_problem_ptr->print_surface_profile(rt_problem_ptr->I_field_, 0, 0, 0, N_theta/2 + 1, 0);     
    rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2 + 1, 0, 1); 
    rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2 + 1, 0, 2);
    rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2 + 1, 0, 3); 


    rt_problem_ptr->print_surface_profile(rt_problem_ptr->I_field_, 0, 2, 2, N_theta/2 + 1, 0);  
    rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 2, 2, N_theta/2 + 1, 0, 1); 
    rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 2, 2, N_theta/2 + 1, 0, 2);
    rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 2, 2, N_theta/2 + 1, 0, 3); 

    // rt_solver.apply_formal_solver();
    // rt_problem_ptr->print_surface_profile(rt_problem_ptr->I_field_, 0, 0, 0, N_theta/2, 0);     

    // rt_solver.compute_emission();
    // rt_problem_ptr->print_profile(rt_problem_ptr->S_field_, 0, 0, 0, 0, N_theta/2, 0);  
  
    // print memory usage 
    const double byte_to_GB = 1.0 / (1000 * 1024 * 1024);

    unsigned long long vm_usage;
    unsigned long long resident_set;
    
    rii::process_mem_usage(vm_usage, resident_set);
    
    if (rt_problem_ptr->mpi_rank_ == 0) std::cout << "Total memory usage (vm_usage) = "     <<  byte_to_GB * vm_usage     << " GB" << std::endl;
    if (rt_problem_ptr->mpi_rank_ == 0) std::cout << "Total memory usage (resident_set) = " <<  byte_to_GB * resident_set << " GB" << std::endl;
    
    rt_problem_ptr->print_PETSc_mem();
        
    if (save_output) rt_problem_ptr->I_field_->write("I_field.raw");          
  }

  Kokkos::finalize();
  PetscFinalize();
  MPI_CHECK(MPI_Finalize());

  return EXIT_SUCCESS;
}


// TODO use Real insted of double

// output
// python ../../sgrid/scripts/transpose_data.py -x 4 -y 4 -z 70 -b 99 -p
// sigma.raw











