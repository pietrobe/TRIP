#include "RT_solver.hpp"


#define NUM_TOL 1e-14

#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define RESET   "\033[0m"

// Helper to get the next "token" regardless of newlines
bool getNextValue(std::ifstream& file, double& value) 
{
    std::string token;
    // Read until a comma or whitespace
    // We use a custom loop to handle the commas manually
    char c;
    while (file.get(c)) {
        if (c == ',' || std::isspace(c)) {
            if (!token.empty()) break; // We found our number
            continue; // Skip leading commas/spaces
        }
        token += c;
    }
    
    if (token.empty()) return false;
    
    try {
        value = std::stod(token);
        return true;
    } catch (...) {
        return false; 
    }
}

bool compare_csv(const std::string& file1, const std::string& file2, double epsilon) 
{
    std::ifstream f1(file1), f2(file2);
    
    if (!f1.is_open()) 
    {
    	std::cerr << "Cannot open " << file1 << "\n";
    	return false;
	}	

	if (!f2.is_open()) 
	{
    	std::cerr << "Cannot open " << file2 << "\n";
    	return false;
	}

    double val1, val2;
    int count = 0;
    double max_diff = 0;

    while (getNextValue(f1, val1) && getNextValue(f2, val2)) 
    {
        count++;
        if (std::fabs(val1 - val2) > max_diff) max_diff = std::fabs(val1 - val2);
        
        if (std::fabs(val1 - val2) > epsilon) 
        {
            std::cout << "Mismatch at value index " << count << "\n";
            std::cout << "File 1: " << val1 << "\n";
            std::cout << "File 2: " << val2 << "\n";
            std::cout << "Absolute difference: " << std::fabs(val1 - val2) << "\n";
            return false;
        }
    }

    // Check if one file has more numbers than the other
    bool more1 = getNextValue(f1, val1);
    bool more2 = getNextValue(f2, val2);

    if (more1 || more2) 
    {
        std::cerr << "Files have different total counts of numbers.\n";
        return false;
    }

    std::cout << "Maximum absolute difference: " << max_diff << std::endl;

    return true;
}

//////////////////////////////////////////////////////////////////////////
// Build executable for testing
//////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
	MPI_CHECK(MPI_Init(&argc, &argv));

	int mpi_rank, mpi_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	
	// import input file config.yml
	const std::string yaml_config_file = getOptionArgument(argc, argv, "--yaml_config");
	AppConfig		  cfg;
	try
	{
		if (mpi_rank == 0) std::cout << "Loading config file: " << yaml_config_file << std::endl;
		cfg = loadConfig(yaml_config_file);
	}
	catch (const std::exception &e)
	{
		std::cerr << "Config error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	PetscInitialize(&argc, &argv, (char *)0, NULL);
	Kokkos::initialize(argc, argv);

#if ACC_SOLAR_3D == _ON_

#pragma message("ACC version ENABLED")

	const int devices_per_node = 4; // for MareNostrum 5 ACC & Daint Alps.
	int		  mpi_error		   = RII_epsilon_contrib::RII_contrib_MPI_Init(devices_per_node, MPI_COMM_WORLD);

	set_RII_contrib_block_size(RII_CONTRIB_BLOCK_SIZE);

	if (RII_epsilon_contrib::RII_contrib_MPI_Is_Device_Handler()) //
		RII_epsilon_contrib::RII_contrib_MPI_Init_Memory_Pool();  //

	const int devices_cnt =									   //
		acc_devices_print_info(mpi_rank, mpi_size, std::cout); //

#else

	set_RII_contrib_block_size(1);

#endif // ACC_SOLAR_3D
	
	{ // start scope for RT_problem and RT_solver		

		// create problem object
		auto rt_problem_ptr = std::make_shared<RT_problem>(cfg);

		// create solver object and solve
		RT_solver rt_solver(rt_problem_ptr);			
    	rt_solver.solve();		
		
		// write output Surface profiles for all surface points
		const int N_x = rt_problem_ptr->N_x_;
		const int N_y = rt_problem_ptr->N_y_;

		const std::string test_out_path  = "../tests/test_output/test_profiles";

		for (int i = 0; i < N_x; ++i)
		{
			for (int j = 0; j < N_y; ++j)
			{		
				rt_problem_ptr->write_surface_point_profiles_csv(test_out_path, i, j);

				if (mpi_rank == 0)
				{
					std::string test_file      = test_out_path + "_" + std::to_string(i) + "_" + std::to_string(j) + ".csv";
					std::string reference_file = (cfg.reference_sol_directory).string() + "profiles_" + std::to_string(i) + "_" + std::to_string(j) + ".csv";        

					if (!compare_csv(test_file, reference_file, NUM_TOL))			
			    	{			    		
			    		std::cout << RED << "==== EROR: TEST NOT PASSED ====\n" << RESET;
			    	}    		
			    	else
			    	{
			    		std::cout << GREEN << "==== TEST PASSED ====\n" << RESET;
			    	}
				}
			}
		}			
	} // end scope for RT_problem and RT_solver

#if ACC_SOLAR_3D == _ON_
	if (RII_epsilon_contrib::RII_contrib_MPI_Is_Device_Handler())
		RII_epsilon_contrib::RII_contrib_MPI_Finalize_Memory_Pool();

	RII_epsilon_contrib::RII_contrib_MPI_Finalize();
#endif // ACC_SOLAR_3D

	Kokkos::finalize();
	PetscFinalize(); 
	MPI_CHECK(MPI_Finalize());

	return 0;
}
