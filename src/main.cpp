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
    const bool save_raw = false;
    const bool use_B    = true;
    const bool use_CRD  = true;
    const bool use_prec = (not use_CRD);

    // inputs
    // const std::string FAL_input_path = "../input/FAL-C/1_B0_V0_12T_8C_99F";
    // const std::string FAL_input_path = "../input/FAL-C/1_B0_V0_12T_8C_64F";
    const std::string FAL_input_path = "../input/FAL-C/96F";

    // const char* PORTA_input_path = "../input/PORTA/cai_0Bx_0By_0Bz_1Vx_1Vy_1Vz_GT4_5x5x133_it100.pmd";
    // const char* PORTA_input_path = "../input/PORTA/cai_1Bx_1By_1Bz_1Vx_1Vy_1Vz_GT4_32x32x133.pmd";
    const char* PORTA_input_path = "../input/PORTA/cai_1Bx_1By_1Bz_1Vx_1Vy_1Vz_GT4_64x64x133.pmd";
  
    auto rt_problem_ptr = std::make_shared<RT_problem>(PORTA_input_path, FAL_input_path, use_CRD, use_B);
    
    const int N_theta = rt_problem_ptr->N_theta_;
    const int N_chi   = rt_problem_ptr->N_chi_; 

    // //FAL-C input    
    // const int N_theta = 8;
    // const int N_chi   = 16;
    // auto rt_problem_ptr = std::make_shared<RT_problem>(FAL_input_path, N_theta, N_chi, use_CRD, use_B);    

    RT_solver rt_solver(rt_problem_ptr, "BESSER", use_prec);
    // RT_solver rt_solver(rt_problem_ptr, "DELO_linear", use_prec);

    rt_solver.solve(); 
   
    // write output
    const std::string output_path = "../output/surface_profiles/";
    const std::string output_file = (use_CRD) ? output_path + "profiles_CRD" : output_path + "profiles_PRD";
    
    const int N_x = rt_problem_ptr->N_x_;
    const int N_y = rt_problem_ptr->N_y_; 

    for (int i = 0; i < N_x; ++i)
    {
      for (int j = 0; j < N_y; ++j)
      {
        rt_problem_ptr->write_surface_point_profiles(output_file, i, j);
      }
    }

    // return EXIT_SUCCESS;
    
    // rt_problem_ptr->print_surface_profile(rt_problem_ptr->I_field_, 0, 0, 0, N_theta/2, 0);     
    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2, 0, 1); 
    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2, 0, 2); 
    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2, 0, 3); 

    // rt_solver.apply_formal_solver();
    // rt_problem_ptr->print_surface_QI_profile(rt_problem_ptr->I_field_, 0, 0, N_theta/2 , 0, 2); 
    // rt_problem_ptr->print_profile(rt_problem_ptr->I_field_, 2, 0, 0, 0, N_theta/2, 0);
    // rt_problem_ptr->print_profile(rt_problem_ptr->I_field_, 2, 0, 0, 0, N_theta/2, 1);    
    // // save_vec(rt_problem_ptr->I_vec_, "../output/I_field_DELO.m" ,"I_DELO");  
    // rt_problem_ptr->I_field_->write("I_field_DELO.raw");          
  
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
        
    if (save_raw) rt_problem_ptr->I_field_->write("I_field.raw");          
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










