#include "RT_solver.hpp"
#include "Test_rii_include.hpp"
#include <chrono>

// compile and run with:
// make -j 32 && srun -n 4 ./main 2 4

int main(int argc, char *argv[]) {

  // return test_solar_3D::test_rii_3D_include(argc, argv);

  MPI_CHECK(MPI_Init(&argc, &argv));
  PetscInitialize(&argc, &argv, (char *)0, NULL);
  Kokkos::initialize(argc, argv);

  {
    const bool using_prec  = false;
    const bool save_output = false;

    const size_t N_theta = atoi(argv[1]);
    const size_t N_chi   = atoi(argv[2]);

    RT_problem rt_problem("../input/FAL-C/1_B0_V0_12T_8C_99F", N_theta, N_chi);

    auto rt_problem_ptr = std::make_shared<RT_problem>(rt_problem);

    RT_solver rt_solver(rt_problem_ptr, "DELO_linear", using_prec);    
    rt_solver.solve();  
    // rt_solver.apply_formal_solver();    
    // rt_solver.compute_emission();
    // rt_solver.test_transfer();

    // const std::string filename =  "../output/pp_" + std::to_string(rt_problem_ptr->mpi_size_) + ".m";
    // const std::string varible  =  "pp" + std::to_string(rt_problem_ptr->mpi_size_);
    // save_vec(rt_problem_ptr->I_vec_, filename.c_str(), varible.c_str());            
    
    rt_problem_ptr->print_surface_profile(rt_problem_ptr->I_field_, 0, 0, 0, N_theta/2, 0);
    rt_problem_ptr->print_surface_profile(rt_problem_ptr->I_field_, 1, 0, 0, N_theta/2, 0);
    // rt_problem_ptr->print_profile(rt_problem_ptr->I_field_, 0, 0, 0, 0, 1, 0);
    
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
