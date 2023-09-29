#include "RT_solver.hpp"
#include "Test_rii_include.hpp"
#include <chrono>
#include "tools.h"

// compile and run with:
// make -j 32 && srun -n 4 ./main

int main(int argc, char *argv[]) {

  MPI_CHECK(MPI_Init(&argc, &argv));
  PetscInitialize(&argc, &argv, (char *)0, NULL);
  Kokkos::initialize(argc, argv);

  {
    const bool save_output = false;
    const bool use_CRD     = false;
    const bool use_prec    = true;
    
    // const std::string input_path  = "../input/PORTA";
    const std::string input_path  = "/home/simone/git/atmos_1d/CaI_4227/tests_problems/B20_V0_12T_8C_99F/";
    // const std::string input_path  = "../input/FAL-C/1_B0_V0_12T_8C_64F_64";
    // const std::string input_path  = "../input/FAL-C/1_B0_V0_12T_8C_4F_64";

    // // PORTA input
    // auto rt_problem_ptr = std::make_shared<RT_problem>("../input/PORTA/cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd", input_path, use_CRD);
    // const size_t N_theta = rt_problem_ptr->N_theta_;
    // const size_t N_chi   = rt_problem_ptr->N_chi_; 

    // FAL-C input
    const size_t N_theta = 4;
    const size_t N_chi   = 4;
    auto rt_problem_ptr = std::make_shared<RT_problem>(input_path, N_theta, N_chi, use_CRD);    

    // RT_solver rt_solver(rt_problem_ptr, "BESSER", use_prec);
    RT_solver rt_solver(rt_problem_ptr, "DELO_linear", use_prec);

    // rt_solver.compute_emission(); ///////////

    // return EXIT_SUCCESS;

    rt_solver.solve(); 

    // rt_problem_ptr->print_surface_profile(rt_problem_ptr->I_field_, 0, 0, 0, N_theta - 1, 0);     
    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta - 1, 0, 1); 

    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2    , 0, 1);
    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2 + 1, 0, 1);
    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta - 1  , 0, 1);
    
    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2    , 0, 2);
    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2 + 1, 0, 2);
    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta - 1  , 0, 2);

    // for (int i = N_theta/2; i < N_theta; ++i)
    // {
      // std::cout << "========== i = " << i << "==========" << std::endl;
      
      rt_problem_ptr->print_surface_profile(rt_problem_ptr->I_field_, 0, 0, 0, N_theta - 1, N_chi - 1);  
      rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta - 1, N_chi - 1, 1);
      rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta - 1, N_chi - 1, 2);
      rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta - 1, N_chi - 1, 3);
    // }
    

    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2    , 0, 3);
    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2 + 1, 0, 3);
    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta - 1  , 0, 3);


    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2 + 1, i_chi, 3); 
   
    // rt_solver.apply_formal_solver(); ////////////
    // // rt_problem_ptr->print_profile(rt_problem_ptr->I_field_, 0, 0, 0, 0, N_theta/2, 0);    
    // // save_vec(rt_problem_ptr->I_vec_, "../output/I_field_DELO.m" ,"I_DELO");  
    // rt_problem_ptr->I_field_->write("I_field_DELO.raw");          
  
    // rt_solver.compute_emission(); ///////////
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










