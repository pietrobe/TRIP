#include <memory>

#include "Test_rii_include.hpp"
#include <rii_emission_coefficient_3D.h>

namespace test_solar_3D {

int test_rii_3D_include(int argc, char *argv[]) {

  MPI_CHECK(MPI_Init(&argc, &argv));
  PetscInitialize(&argc, &argv, (char *)0, NULL);
//   Kokkos::initialize(argc, argv);

  {
    const int N_theta = atoi(argv[1]);
    const int N_chi = atoi(argv[2]);

    RT_problem rt_problem("/users/sriva/git_rii_tests/r_iii_exact_tests/"
                          "CaI_4227/FAL-C/B0_V0_12T_8C_99F/",
                          N_theta, N_chi);

    auto rt_problem_ptr = std::make_shared<RT_problem>(rt_problem);

    auto rt_solver_sh_ptr =
        std::make_shared<RT_solver>(rt_problem_ptr, "DELO_linear");

    using rii_eps_comp_3D = rii_include::emission_coefficient_computation_3D;
    using rii_formal_solver_factory =
        rii_include::formal_solver_factory_from_3D_RT_problem;

    auto ecc_sh_ptr =
        rii_eps_comp_3D::make_emission_coefficient_computation_3D_shared_ptr();

    auto fsf_sh_ptr = rii_formal_solver_factory::
        make_formal_solver_factory_from_3D_RT_problem_shared_ptr();

    using in_RT_problem_3D =
        rii_include::in_RT_problem_3D_interface<RT_problem>;

    std::cout << "Start add models 3D\n";
    std::cout.flush();

    in_RT_problem_3D::add_models(rt_problem_ptr, ecc_sh_ptr, fsf_sh_ptr, true);
    std::cout << "End add models 3D\n";
    std::cout.flush();

    std::cout << "Start make FS 3D\n";
    std::cout.flush();

    fsf_sh_ptr->make_formal_solver();

    std::cout << "End make FS 3D\n";
    std::cout.flush();

    using emission_coefficient_components = rii_include::
        emission_coefficient_computation::emission_coefficient_components;

    std::cout << "Start make_computation_function\n";
    std::cout.flush();

    std::list<emission_coefficient_components> components{
        emission_coefficient_components::epsilon_R_II,
        emission_coefficient_components::epsilon_R_III_GL,
        emission_coefficient_components::epsilon_csc};

    std::cout << fsf_sh_ptr->get_formal_solver_shared_ptr()
                     ->print_emission_model(components);

    const auto epsilon_fun = ecc_sh_ptr->make_computation_function(components);

    const auto offset_fun = rii_include::make_default_offset_function(
        N_theta, N_chi,
        ecc_sh_ptr->get_formal_solver()
            ->get_solar_atmosphere_model()
            ->get_frequencies_grid_size());

    std::vector<double> input_array(90000);

    const auto in_field = ecc_sh_ptr->update_incoming_field(2, 2, 3, offset_fun,
                                                            input_array.data());

    const auto out_field = epsilon_fun(2, 2, 3);

    if (not out_field) {
      std::cout << "if (not out_field) { \n";
      exit(0);
    }

    std::vector<double> out_array(90000);
    rii_include::make_indices_convertion_function<double>(
        out_field, offset_fun)(out_array.data());

    std::cout << "End make_computation_function\n";
    std::cout.flush();

    // rt_solver_sh_ptr->solve();
  }

//   Kokkos::finalize();
  PetscFinalize();
  MPI_CHECK(MPI_Finalize());

  return EXIT_SUCCESS;
}

} // namespace test_solar_3D