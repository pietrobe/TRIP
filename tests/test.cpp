#ifndef TEST_MAIN_CPP
#define TEST_MAIN_CPP

#include "RT_solver.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>

#define ABS_TOL 1e-7
#define REL_TOL 1e-6

#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define RESET   "\033[0m"


// Trim leading and trailing whitespace
std::string trim(const std::string& str)
{
    size_t start = str.find_first_not_of(" \t\r\n");
    size_t end = str.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}


// Split CSV line into cells, trimming whitespace
std::vector<std::string> split_csv_line(const std::string& line)
{
    std::vector<std::string> cells;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) {
        cells.push_back(trim(cell));
    }
    return cells;
}

// Compare two CSV files
std::string compare_csv_report(const std::string& file1, const std::string& file2, bool &is_equal)
{
    std::ifstream f1(file1), f2(file2);
    if (!f1.is_open()) { is_equal = false; return "Cannot open " + file1; }
    if (!f2.is_open()) { is_equal = false; return "Cannot open " + file2; }

    const char* labels[4] = {"I","Q","U","V"};
    const double absTol[4] = {ABS_TOL, ABS_TOL, ABS_TOL, ABS_TOL};
    const double relTol[4] = {REL_TOL, REL_TOL*100, REL_TOL*100, REL_TOL*100};

    double max_abs_diff[4] = {0.0,0.0,0.0,0.0};
    double max_rel_diff[4] = {0.0,0.0,0.0,0.0};

    std::string line1, line2;
    int row = 0;
    bool anyMismatch = false;
    const double small_value_eps = 1e-8; 

    std::ostringstream report;
    report << "Comparing files:\n"
           << "  File1: " << file1 << "\n"
           << "  File2: " << file2 << "\n\n";

    while (std::getline(f1, line1) && std::getline(f2, line2))
    {
        row++;
        int stokesIndex = (row-1) % 4;
        double abs_epsilon = absTol[stokesIndex];
        double rel_epsilon = relTol[stokesIndex];

        auto cells1 = split_csv_line(line1);
        auto cells2 = split_csv_line(line2);

        if (cells1.size() != cells2.size()) {
            report << "Row " << row << ": Column count mismatch (" 
                   << cells1.size() << " vs " << cells2.size() << ")\n\n";
            anyMismatch = true;
            continue;
        }

        bool row_has_mismatch = false;
        std::ostringstream row_report;

        for (size_t col = 0; col < cells1.size(); ++col)
        {
            double val1 = std::stod(cells1[col]);
            double val2 = std::stod(cells2[col]);

            double abs_error = std::fabs(val1 - val2);
            double denom = std::max(std::fabs(val1), std::fabs(val2));
            double rel_error = (denom == 0.0) ? 0.0 : abs_error / denom;

            bool small_values_check = (std::fabs(val1) > small_value_eps || std::fabs(val2) > small_value_eps);
            bool mismatch = (abs_error > abs_epsilon) || (small_values_check && rel_error > rel_epsilon);

            max_abs_diff[stokesIndex] = std::max(max_abs_diff[stokesIndex], abs_error);
            max_rel_diff[stokesIndex] = std::max(max_rel_diff[stokesIndex], rel_error);

            if (mismatch) {
                row_has_mismatch = true;
                row_report << "    Col " << col+1 << ": "
                           << "File1=" << val1 << ", "
                           << "File2=" << val2 << ", "
                           << "AbsErr=" << abs_error << ", "
                           << "RelErr=" << rel_error << "\n";
            }
        }

        if (row_has_mismatch) {
            anyMismatch = true;
            report << "Row " << row << " (Variable " << labels[stokesIndex] << "):\n";
            report << row_report.str() << "\n";
        }
    }

    if (std::getline(f1, line1) || std::getline(f2, line2)) {
        report << "Row count mismatch detected\n\n";
        anyMismatch = true;
    }

    report << "** Maximum Differences per Variable **\n";
    report << "Variable | MaxAbsDiff      | MaxRelDiff\n";
    report << "---------|----------------|----------------\n";
    for (int i=0; i<4; i++)
    {
        report << labels[i] << "        | "
               << max_abs_diff[i] << " | "
               << max_rel_diff[i] << "\n";
    }
    report << "\n";

    is_equal = !anyMismatch;
    return report.str();
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


                    bool equal;
                    std::string report = compare_csv_report(test_file, reference_file, equal);
                    if (!equal) {
                        std::cout << RED << "==== ERROR: TEST NOT PASSED ====\n" << RESET;
                    } else {
                        std::cout << GREEN << "==== TEST PASSED ====\n" << RESET;
                    }
                    std::cout << report;
				}
			}
		}			
	} // end scope for RT_problem and RT_solver

#if ACC_SOLAR_3D == _ON_
	if (RII_epsilon_contrib::RII_contrib_MPI_Is_Device_Handler())
		RII_epsilon_contrib::RII_contrib_MPI_Finalize_Memory_Pool();

	RII_epsilon_contrib::RII_contrib_MPI_Finalize();
#endif // ACC_SOLAR_3D

	PetscFinalize(); 
	MPI_CHECK(MPI_Finalize());

	return 0;
}

#endif // TEST_MAIN_CPP