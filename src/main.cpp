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
    const bool output   = false;
    const bool use_B    = false;
    const bool use_CRD  = false;
    const bool use_prec = (not use_CRD);

    // inputs
    const std::string FAL_input_path = "/users/sriva/git/solar_3d/input/FAL-C/B20_V0_12T_8C_99F_1Pi4_9Pi8";
    // const std::string FAL_input_path = "../input/FAL-C/1_B0_V0_12T_8C_64F";
    // const std::string FAL_input_path = "/users/pietrob/solar_3d/input/FAL-C/96F";
    // const std::string FAL_input_path = "/users/pietrob/solar_3d/input/FAL-C/64F";

    // const char* PORTA_input_path = "/users/pietrob/solar_3d/input/PORTA/cai_0Bx_0By_0Bz_1Vx_1Vy_1Vz_GT4_5x5x133_it100.pmd";
    // // const char* PORTA_input_path = "/users/pietrob/solar_3d/input/PORTA/cai_1Bx_1By_1Bz_1Vx_1Vy_1Vz_GT4_32x32x133.pmd";
    // // const char* PORTA_input_path = "/users/pietrob/solar_3d/input/PORTA/cai_1Bx_1By_1Bz_1Vx_1Vy_1Vz_GT4_64x64x133.pmd";
  
    // auto rt_problem_ptr = std::make_shared<RT_problem>(PORTA_input_path, FAL_input_path, use_CRD, use_B);
    
    // const int N_theta = rt_problem_ptr->N_theta_;
    // const int N_chi   = rt_problem_ptr->N_chi_; 

    //FAL-C input    
    const int N_theta = 8;
    const int N_chi   = 16;
    auto rt_problem_ptr = std::make_shared<RT_problem>(FAL_input_path, N_theta, N_chi, use_CRD, use_B);    

    RT_solver rt_solver(rt_problem_ptr, "BESSER", use_prec);
    // RT_solver rt_solver(rt_problem_ptr, "DELO_linear", use_prec);

    rt_solver.solve(); 
    // rt_solver.solve_checkpoint("../output/surface_profiles_5x5x133/", 20); 
    
    // write output
    if (output)
    {            
        const std::string output_path = "../output/surface_profiles_test/"; // TODO change
        // const std::string output_path = "output/surface_profiles_64x64x133/B0V0/"; // TODO change
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

        // free some memory    
        rt_problem_ptr->free_fields_memory(); 
        rt_solver.free_fields_memory();

        // write is arbitriary direction     
        Real chi   = 0.19635;
        Real mu    = 0.930568;
        Real theta = acos(mu);

        std::string output_file_Omega;

        rt_solver.apply_formal_solver_Omega(theta, chi);
        
        output_file_Omega = output_file + "_control";

        // for (int i = 0; i < N_x; ++i)
        // {
        //    for (int j = 0; j < N_y; ++j)
        //    {
              rt_problem_ptr->write_surface_point_profiles_Omega(output_file_Omega, 0, 0);
        //    }
        // }

        mu = 0.1;
        theta = acos(mu);

        // allocate new data structure and compute I_Field_Omega
        rt_solver.apply_formal_solver_Omega(theta, chi);

        output_file_Omega = output_file + "_mu01";
        for (int i = 0; i < N_x; ++i)
        {
           for (int j = 0; j < N_y; ++j)
           {
                rt_problem_ptr->write_surface_point_profiles_Omega(output_file_Omega, i, j);
           }
        }

        mu = 1.0;
        theta = acos(mu);

        // allocate new data structure and compute I_Field_Omega
        rt_solver.apply_formal_solver_Omega(theta, chi);

        output_file_Omega = output_file + "_mu1";

        for (int i = 0; i < N_x; ++i)
        {
           for (int j = 0; j < N_y; ++j)
           {
                rt_problem_ptr->write_surface_point_profiles_Omega(output_file_Omega, i, j);
           }
        }

        // if (save_raw) rt_problem_ptr->I_field_->write("/scratch/snx3000/pietrob/I_field.raw");          
          
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
    }
  
    // print memory usage 
    const double byte_to_GB = 1.0 / (1000 * 1024 * 1024);

    unsigned long long vm_usage;
    unsigned long long resident_set;
    
    rii::process_mem_usage(vm_usage, resident_set);
    
    if (rt_problem_ptr->mpi_rank_ == 0) std::cout << "Total memory usage (vm_usage) = "     <<  byte_to_GB * vm_usage     << " GB" << std::endl;
    if (rt_problem_ptr->mpi_rank_ == 0) std::cout << "Total memory usage (resident_set) = " <<  byte_to_GB * resident_set << " GB" << std::endl;
    
    rt_problem_ptr->print_PETSc_mem();        
  }
  
  Kokkos::finalize();
  PetscFinalize(); //CHKERRQ(ierr);
  MPI_CHECK(MPI_Finalize());

  return 0;
}


// TODO use Real insted of double

// output
// python ../../sgrid/scripts/transpose_data.py -x 4 -y 4 -z 70 -b 99 -p
// sigma.raw










