#include "RT_solver.hpp"
#include "Test_rii_include.hpp"
#include <chrono>

//  make -j4 && mpirun -n 4 ./main

int main(int argc, char *argv[]) {

  // return test_solar_3D::test_rii_3D_include(argc, argv);
  

  MPI_CHECK(MPI_Init(&argc, &argv));
  PetscInitialize(&argc, &argv, (char *)0, NULL);
  Kokkos::initialize(argc, argv);

  {
    const size_t N_theta = atoi(argv[1]);
    const size_t N_chi = atoi(argv[2]);

    RT_problem rt_problem("../input/FAL-C/1_B0_V0_12T_8C_99F", N_theta, N_chi);

    auto rt_problem_ptr = std::make_shared<RT_problem>(rt_problem);

    RT_solver rt_solver(rt_problem_ptr, "DELO_linear");

    rt_solver.solve();
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
