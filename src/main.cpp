#include "RT_solver.hpp"
#include "Test_rii_include.hpp"
#include "RT_utility.hpp"

#include <chrono>
#include "tools.h"
#include <string>
#include <filesystem>
#include <iomanip>
#include <sstream>

// WARNING: if you want to use PORTA input for 3D setup, you need to set USE_PORTA_INPUT = 1
// otherwise, it will use FAL-C input for 1D setup
#define USE_PORTA_INPUT 1
//////////////////////////////////////////////////////////////////////////

// WARNING: if you want to use command line options, you need to set USE_CMD_LINE_OPTIONS = 1
// otherwise, it will use the default and hard-coded values
#define USE_CMD_LINE_OPTIONS 1
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {

  std::stringstream ss_a, ss_b;
  std::filesystem::path output_info_file;

  MPI_CHECK(MPI_Init(&argc, &argv));

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  { // check if the user wants to print the help message
    // if yes, print the help message and return
    
    if (getOptionFlag(argc, argv, "--help")) {
       if (mpi_rank == 0) 
          print_help();
          
      return 0;
    }
  }

  PetscInitialize(&argc, &argv, (char *)0, NULL);
  Kokkos::initialize(argc, argv);

  {
    const bool output   = true;
    const bool output_overwrite_prevention = true; // if true the application stops (with an MPI_Abort) if the output directory already exists

    const bool use_B    = true;

#if USE_CMD_LINE_OPTIONS == 1
    const bool use_CRD  = getOptionFlag(argc, argv, "--CRD");
    const bool use_prec = (not use_CRD);
#else
    const bool use_CRD  = true;
    const bool use_prec = (not use_CRD);
#endif


  // Set here the main input and output directories //////////////////////////
#if USE_CMD_LINE_OPTIONS == 1
  const std::filesystem::path main_input_dir  = getOptionArgument(argc, argv, "--input_dir");
  const std::filesystem::path main_output_dir = getOptionArgument(argc, argv, "--output_dir");
#else
//  const std::filesystem::path main_input_dir  = "../input/PORTA";
  const std::filesystem::path main_input_dir  = "/users/pietrob/solar_3d/input/PORTA";  
  const std::filesystem::path main_output_dir = "/users/pietrob/solar_3d/output";
#endif
  ////////////////////////////////////////////////////////////////////////////

#if USE_PORTA_INPUT == 1   // PORTA setup for 3D

  // Set here the problem input file //////////////////////////
  #if USE_CMD_LINE_OPTIONS == 1
    std::string input_pmd_string = getOptionArgument(argc, argv, "--problem_pmd_file");
    
    std::string input_cul_string;
    std::string input_qel_string;
    std::string input_llp_string;
    std::string input_back_string;

    const std::string input_config_string = getOptionArgument(argc, argv, "--problem_input_config");

    if (not input_config_string.empty()) {
      
      const auto config_map = get_input_PORTA(main_input_dir / std::filesystem::path(input_config_string), mpi_rank);

      input_pmd_string = config_map.at("pmd");

      input_cul_string = config_map.at("cul");
      input_qel_string = config_map.at("qel");
      input_llp_string = config_map.at("llp");
      input_back_string = config_map.at("back");
    }

  #else
    const std::string input_pmd_string = std::string("AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_V0_conv.pmd");
    const std::string input_llp_string = std::string("AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_V0_conv.llp");
    
    const std::string input_cul_string  = std::string("AR_385_Cut_64x64_mirrorxy-CRD_I_V0.cul");
    const std::string input_qel_string  = std::string("AR_385_Cut_64x64_mirrorxy-CRD_I_V0.qel");
    const std::string input_back_string = std::string("AR_385_Cut_64x64_mirrorxy-CRD_I_V0.back");
  #endif
  /////////////////////////////////////////////////////////////

    // const auto input_pmd_file = std::filesystem::path(input_pmd_string);

    auto frequencies_input_path =  main_input_dir / std::filesystem::path("frequency/96F");

    auto PORTA_input_pmd =  main_input_dir / std::filesystem::path(input_pmd_string);

    // lambda to build the RT_problem object
    auto create_rt_problem = [&]() {

        if (input_cul_string.empty() or input_qel_string.empty() or input_llp_string.empty()) {
        // VEECHIO  // solo PMD input at least one of the cul, qel, llp is missing
            if (mpi_rank == 0) {
              std::cout << "WARNING: using ONLY PMD input file" << std::endl;
            }

            return std::make_shared<RT_problem>(PORTA_input_pmd.string().c_str(), frequencies_input_path.string(), use_CRD, use_B);

        } else {

          if (mpi_rank == 0) {
            std::cout << "WARNING: using PMD + CUL + QEL + LLP + BACK input files" << std::endl;
          }

        // NUOVO // PMD + CUL + QEL + LLP input
            // create cul and qel input path
            auto input_cul_path  = main_input_dir / std::filesystem::path(input_cul_string);
            auto input_qel_path  = main_input_dir / std::filesystem::path(input_qel_string);
            auto input_llp_path  = main_input_dir / std::filesystem::path(input_llp_string);
            auto input_back_path = main_input_dir / std::filesystem::path(input_back_string);

            return std::make_shared<RT_problem>(PORTA_input_pmd.string().c_str(),
                                                input_cul_path.string().c_str(),
                                                input_qel_path.string().c_str(),
                                                input_llp_path.string().c_str(),
                                                input_back_path.string().c_str(),
                                                frequencies_input_path.string(), use_CRD, use_B);
        }
    }; // end lambda create_rt_problem

    auto rt_problem_ptr = create_rt_problem();  /// CALL LAMBDA to create RT_problem object

    const int N_theta = rt_problem_ptr->N_theta_;
    const int N_chi   = rt_problem_ptr->N_chi_; 

    if (rt_problem_ptr->mpi_rank_ == 0) {

      // print the command line arguments
      ss_a << "MPI size = " << mpi_size << std::endl;
      ss_a << "Date and time: " << getCurrentDateTime() << std::endl; 
      ss_a << std::endl << "Command line arguments: " << std::endl;
      ss_a << "argc: " << argc << std::endl;
      ss_a << "argv: ";
      for (int i = 0; i < argc; ++i) {
        ss_a << argv[i] << " ";
      }
      ss_a << std::endl << std::endl;

      ss_a << "PORTA 3D input file: " << PORTA_input_pmd << std::endl;
      ss_a << "Frequencies input path: " << frequencies_input_path << std::endl;
      
      ss_a << "input_cul_string: " << (input_cul_string.empty() ? std::string("not provided") : input_cul_string) << std::endl;
      ss_a << "input_qel_string: " << (input_qel_string.empty() ? std::string("not provided") : input_qel_string) << std::endl;
      ss_a << "input_llp_string: " << (input_llp_string.empty() ? std::string("not provided") : input_llp_string) << std::endl;

      ss_a << "N_theta =            " << N_theta << std::endl;
      ss_a << "N_chi =              " << N_chi << std::endl << std::endl;

      const auto petsc_index_size = sizeof(PetscInt);
      ss_a << "PetscInt size: " << petsc_index_size << " bytes; " << (petsc_index_size * 8) << " bits." << std::endl << std::endl;

      std::cout << ss_a.str();  
    }

    

#else 
    //FAL-C input for 1D input setup

    // inputs
    auto problem_input_file = std::filesystem::path("FAL-C/B20_V0_12T_8C_99F_1Pi4_9Pi8");
    auto problem_input_FAL = main_input_dir / problem_input_file;
    // const std::string FAL_input_path = "../input/FAL-C/1_B0_V0_12T_8C_64F";
    // const std::string FAL_input_path = "/users/pietrob/solar_3d/input/FAL-C/96F";
    // const std::string FAL_input_path = "/users/pietrob/solar_3d/input/FAL-C/64F";

    const int N_theta = 8;
    const int N_chi   = 16;
    auto rt_problem_ptr = std::make_shared<RT_problem>(problem_input_FAL.string(), N_theta, N_chi, use_CRD, use_B);    
#endif

    RT_solver rt_solver(rt_problem_ptr, "BESSER", use_prec);
   // RT_solver rt_solver(rt_problem_ptr, "DELO_linear", use_prec);

    //////////////////////////////////////////////////////////////////////////
    // Prepare output directory
    // If the output directory does not exist, create it
    // It it exists, abort !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    std::string output_file;
    if (output)
    {            
        const std::filesystem::path output_path = main_output_dir / std::filesystem::path(input_pmd_string + ((use_CRD) ? ".CRD" : ".PRD"));

        //  if (rt_problem_ptr->mpi_rank_ == 0) 
        //    std::cout << "Output path: " << output_path << std::endl;

        if (rt_problem_ptr->mpi_rank_ == 0) {
          if (not std::filesystem::exists(output_path)){
            std::filesystem::create_directories(output_path);
          } else if (output_overwrite_prevention) {
            std::cerr << "File: " << __FILE__ << " Line: " << __LINE__ << " MPI rank: " << rt_problem_ptr->mpi_rank_ << " Directory: " << output_path << " already exists" << std::endl;
            std::cerr << "Use different output directory" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
          }
        } 

        // const std::string output_path = "output/surface_profiles_64x64x133/B0V0/"; // TODO change
        output_file = (use_CRD) ? (output_path / "profiles_CRD").string() : (output_path / "profiles_PRD").string();

       if (rt_problem_ptr->mpi_rank_ == 0)  {
        ss_b << "Output directory: " << output_path << std::endl;
        std::cout << ss_b.str();

        const auto output_info_file = output_path / "info.txt";
        std::ofstream output_file_info(output_info_file);

        output_file_info << ss_a.str();
        output_file_info << ss_b.str();
        output_file_info.close();
      }
    }

    ///////////////////////////////////////////////////
    // solve //////////////////////////////////////////
    rt_solver.solve(); 
    // rt_solver.solve_checkpoint("../output/surface_profiles_5x5x133/", 20); 
    
    // lambda to compute arbitrary beam
      const auto compute_arbitrary_beam = [&] (const Real mu, const Real chi, const std::string output_file) {

        std::string output_file_Omega_mu = output_file + "_mu" + std::to_string(mu);
        
        const Real theta = acos(mu);
        rt_solver.apply_formal_solver_Omega(theta, chi);

        for (int i = 0; i < N_theta; ++i)
        {
          for (int j = 0; j < N_chi; ++j)
          {
                rt_problem_ptr->write_surface_point_profiles_Omega(output_file_Omega_mu, i, j);
          }
        }
    };


    // write output
    if (output){
        const int N_x = rt_problem_ptr->N_x_;
        const int N_y = rt_problem_ptr->N_y_; 

        for (int i = 0; i < N_x; ++i)
        {
           for (int j = 0; j < N_y; ++j)
           {
            rt_problem_ptr->write_surface_point_profiles(output_file, i, j);
           }
        }      

        // rt_problem_ptr->write_surface_point_profiles(output_file, 0, 0);
        
        // old code: copied below .....

        // // std::vector<Real> mus = {0.1, 1.0}; //// ATTENTION: arbitrary beam directions
        std::vector<Real> mus = {}; //// ATTENTION: arbitrary beam directions NO ARBITRARY BEAMS
        Real chi   = 0.19635;
    
        if (rt_problem_ptr->mpi_rank_ == 0 and mus.size() ==0 ){
          std::cout << "WARNING: no arbitrary beams" << std::endl;
        } else if (rt_problem_ptr->mpi_rank_ == 0) {
          std::cout << "Arbitrary beams: ";
          for (auto mu : mus) {
            std::cout << mu << " ";
          }
          std::cout << std::endl;
        }
         
        for (Real mu : mus)
        {
          compute_arbitrary_beam(mu, chi, output_file);
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
    } // end write output
  
    // print memory usage 
    const double byte_to_GB = 1.0 / (1000 * 1024 * 1024);

    unsigned long long vm_usage;
    unsigned long long resident_set;
    
    rii::process_mem_usage(vm_usage, resident_set);
    
    if (rt_problem_ptr->mpi_rank_ == 0){
      std::stringstream ss_mem;
      ss_mem << "Total memory usage (vm_usage) = "     <<  byte_to_GB * vm_usage     << " GB" << std::endl;
      ss_mem << "Total memory usage (resident_set) = " <<  byte_to_GB * resident_set << " GB" << std::endl;
    
      std::string mem_petsc = rt_problem_ptr->print_PETSc_mem();        
      ss_mem << mem_petsc << std::endl;

      std::cout << ss_mem.str();

      
      std::ofstream output_file_info(output_info_file, std::ios::app);
      output_file_info << ss_mem.str();
      output_file_info.close();
      
    }
  }
  
  Kokkos::finalize();
  PetscFinalize(); //CHKERRQ(ierr);
  MPI_CHECK(MPI_Finalize());

  return 0;
}


// TODO use Real instead of double

// output
// python ../../sgrid/scripts/transpose_data.py -x 4 -y 4 -z 70 -b 99 -


/// il vecchio codice per scrivere i profili 
/// e per calcolare i profili in una direzione arbitraria

        // // free some memory    
        // rt_problem_ptr->free_fields_memory(); 
        // rt_solver.free_fields_memory();

        // std::string output_file_Omega;

        // rt_solver.apply_formal_solver_Omega(theta, chi);
        
        // output_file_Omega = output_file + "_control";

        // // for (int i = 0; i < N_x; ++i)
        // // {
        // //    for (int j = 0; j < N_y; ++j)
        // //    {
        //       rt_problem_ptr->write_surface_point_profiles_Omega(output_file_Omega, 0, 0);
        // //    }
        // // }

        // mu = 0.1;
        // theta = acos(mu);

        // // allocate new data structure and compute I_Field_Omega
        // rt_solver.apply_formal_solver_Omega(theta, chi);

        // output_file_Omega = output_file + "_mu01";
        // for (int i = 0; i < N_x; ++i)
        // {
        //    for (int j = 0; j < N_y; ++j)
        //    {
        //         rt_problem_ptr->write_surface_point_profiles_Omega(output_file_Omega, i, j);
        //    }
        // }

        // mu = 1.0;
        // theta = acos(mu);

        // // allocate new data structure and compute I_Field_Omega
        // rt_solver.apply_formal_solver_Omega(theta, chi);

        // output_file_Omega = output_file + "_mu1";

        // for (int i = 0; i < N_x; ++i)
        // {
        //    for (int j = 0; j < N_y; ++j)
        //    {
        //         rt_problem_ptr->write_surface_point_profiles_Omega(output_file_Omega, i, j);
        //    }
        // }
