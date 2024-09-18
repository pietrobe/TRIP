#include "RT_solver.hpp"
#include "Test_rii_include.hpp"
#include <chrono>
#include "tools.h"
#include <string>
#include <filesystem>
#include <iomanip>
#include <sstream>

std::string getCurrentDateTime() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

// compile and run with:
// make -j 32 && srun -n 4 ./main

// WARNING: if you want to use PORTA input for 3D setup, you need to set USE_PORTA_INPUT = 1
// otherwise, it will use FAL-C input for 1D setup
#define USE_PORTA_INPUT 1
//////////////////////////////////////////////////////////////////////////

// WARNING: if you want to use command line options, you need to set USE_CMD_LINE_OPTIONS = 1
// otherwise, it will use the default and hard-coded values
#define USE_CMD_LINE_OPTIONS 1
//////////////////////////////////////////////////////////////////////////

std::string getOptionArgument(int argc, char *argv[], const std::string &option) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == option && i + 1 < argc) {
      return argv[i + 1];
    }
  }
  return std::string();  // Option not found or argument missing
}

bool getOptionFlag(int argc, char *argv[], const std::string &option) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == option) {
      return true;
    }
  }
  return false;
}

void print_help() {
  std::cout << "----------------------------------------------------------------" << std::endl << std::endl;
  std::cout << "Usage: mpirun ./solar_3D [options] [PETSC options]" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  --CRD: use CRD" << std::endl;
  std::cout << "  --input_dir <input_dir>: input directory" << std::endl;
  std::cout << "  --output_dir <output_dir>: output directory" << std::endl;
  std::cout << "  --problem_input_file <problem_input_file>: problem input file" << std::endl;
  std::cout << "  --help: print help and exit" << std::endl << std::endl;

  std::cout << "----------------------------------------------------------------" << std::endl << std::endl;

  std::cout << "Example: mpirun ./solar_3D --CRD --input_dir /path/to/input --output_dir /path/to/output --problem_input_file problem_input_file.pmd  -ksp_type gmres -ksp_max_it 100 -ksp_rtol 1e-10" << std::endl;
  std::cout << std::endl;
  std::cout << "In the output directory, the code creates a results directory with the name of the problem input file and the extension \".CRD\" or \".PRD\", depending on the --CRD option." << std::endl;
  std::cout << "If the results output directory already exists, the code will stop." << std::endl;
  std::cout << "Default solver is the PRD." << std::endl;
  std::cout << std::endl;
}




//////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {

  std::stringstream ss_a, ss_b;

  MPI_CHECK(MPI_Init(&argc, &argv));

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  { // check if the user wants to print the help message
    // if yes, print the help message and return
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
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
  const std::filesystem::path main_input_dir = getOptionArgument(argc, argv, "--input_dir");
  const std::filesystem::path  main_output_dir = getOptionArgument(argc, argv, "--output_dir");
#else
  const std::filesystem::path main_input_dir = "../input/PORTA";
  const std::filesystem::path main_output_dir = "output";
#endif
  ////////////////////////////////////////////////////////////////////////////

#if USE_PORTA_INPUT == 1   // PORTA setup for 3D

  // Set here the problem input file //////////////////////////
  #if USE_CMD_LINE_OPTIONS == 1
    const std::string problem_input_file_string = getOptionArgument(argc, argv, "--problem_input_file");
  #else
    const std::string problem_input_file_string = std::string("cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd");
  #endif
  /////////////////////////////////////////////////////////////

//    auto PORTA_input_path = main_input_dir / std::filesystem::path("PORTA/cai_0Bx_0By_0Bz_1Vx_1Vy_1Vz_GT4_5x5x133_it100.pmd");
    // auto PORTA_input_path =  main_input_dir / std::filesystem::path("PORTA/cai_1Bx_1By_1Bz_1Vx_1Vy_1Vz_GT4_32x32x133.pmd");
    // auto PORTA_input_path =  main_input_dir / std::filesystem::path("PORTA/cai_1Bx_1By_1Bz_1Vx_1Vy_1Vz_GT4_64x64x133.pmd");

    const auto problem_input_file = std::filesystem::path(problem_input_file_string);

    auto frequencies_input_path =  main_input_dir / std::filesystem::path("frequency/96F");

//    auto PORTA_input_path =  main_input_dir / std::filesystem::path("PORTA/cai_0Bx_0By_0Bz_1Vx_1Vy_1Vz_GT4_5x5x133_it100.pmd");
    // // const char* PORTA_input_path = "/users/pietrob/solar_3d/input/PORTA/cai_1Bx_1By_1Bz_1Vx_1Vy_1Vz_GT4_32x32x133.pmd";
    // auto PORTA_input_path =  main_input_dir / std::filesystem::path("PORTA/cai_1Bx_1By_1Bz_1Vx_1Vy_1Vz_GT4_64x64x133.pmd");
    // auto PORTA_input_path =  main_input_dir / std::filesystem::path("PORTA/cai_1Bx_1By_1Bz_1Vx_1Vy_1Vz_GT4_32x32x133.pmd");
    // auto PORTA_input_path =  main_input_dir / std::filesystem::path("PORTA/cai_1Bx_1By_1Bz_1Vx_1Vy_1Vz_GT4_16x16x133.pmd");

    auto PORTA_input_path =  main_input_dir / problem_input_file;

    auto rt_problem_ptr = std::make_shared<RT_problem>(PORTA_input_path.string().c_str(), frequencies_input_path.string() , use_CRD, use_B);
    
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

      ss_a << "PORTA 3D input file: " << PORTA_input_path << std::endl;
      ss_a << "N_theta =            " << N_theta << std::endl;
      ss_a << "N_chi =              " << N_chi << std::endl << std::endl;

      const auto petsc_index_size = sizeof(PetscInt);
      ss_a << "PetscInt size: " << petsc_index_size << " bytes; " << (petsc_index_size * 8) << " bits." << std::endl << std::endl;

      std::cout << ss_a.str();  
    }

    

#else 
    //FAL-C input for 1D setup

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
        const std::filesystem::path output_path = main_output_dir / std::filesystem::path(problem_input_file_string + ((use_CRD) ? ".CRD" : ".PRD"));

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
// python ../../sgrid/scripts/transpose_data.py -x 4 -y 4 -z 70 -b 99 -
