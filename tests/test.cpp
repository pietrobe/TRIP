#ifndef TEST_MAIN_CPP
#define TEST_MAIN_CPP

#include "RT_solver.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip> 

#define ABS_TOL 1e-5
#define REL_TOL 1e-6

#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define RESET   "\033[0m"

struct CSVMaxDiff {
    double max_abs[4] = {0,0,0,0};
    double max_rel[4] = {0,0,0,0};
};

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
CSVMaxDiff compare_csv_report(const std::string& file1, const std::string& file2, bool &is_equal, bool verbose)
{
    std::ifstream f1(file1), f2(file2);
    if (!f1.is_open()) { 
        is_equal = false; 
        std::cout << "Cannot open " << file1 << "\n"; 
        return CSVMaxDiff{}; 
    }
    if (!f2.is_open()) { 
        is_equal = false; 
        std::cout << "Cannot open " << file2<< "\n";
        return CSVMaxDiff{}; 
    }


    const char*  labels[4] = {"I","Q/I[%]","U/I[%]","V/I[%]"};
    const double absTol[4] = {ABS_TOL, ABS_TOL, ABS_TOL, ABS_TOL};
    const double relTol[4] = {REL_TOL, REL_TOL*100, REL_TOL*100, REL_TOL*100};

    double max_abs_diff[4] = {0.0,0.0,0.0,0.0};
    double max_rel_diff[4] = {0.0,0.0,0.0,0.0};

    std::string line1, line2;
    int row = 0;
    bool anyMismatch = false;
    const double small_value_eps = 1e-8; 

    std::ostringstream report;
    if (verbose) {
        report << "Comparing files:\n"
            << "  File1: " << file1 << "\n"
            << "  File2: " << file2 << "\n\n";
    }

    while (std::getline(f1, line1) && std::getline(f2, line2))
    {
        row++;
        int stokesIndex = (row-1) % 4;
        double abs_epsilon = absTol[stokesIndex];
        double rel_epsilon = relTol[stokesIndex];

        auto cells1 = split_csv_line(line1);
        auto cells2 = split_csv_line(line2);

        if (cells1.size() != cells2.size()) {
            if (verbose) {
                report << "Row " << row << ": Column count mismatch (" 
                    << cells1.size() << " vs " << cells2.size() << ")\n\n";
            }
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
                if (verbose){
                    row_report << "    Col " << col+1 << ": "
                            << "File1=" << val1 << ", "
                            << "File2=" << val2 << ", "
                            << "AbsErr=" << abs_error << ", "
                            << "RelErr=" << rel_error << "\n";
                }
            }
        }

        if (row_has_mismatch) {
            anyMismatch = true;
            if (verbose) {
                report << "Row " << row << " (Variable " << labels[stokesIndex] << "):\n";
                report << row_report.str() << "\n";
            }
        }
    }

    CSVMaxDiff results;
    for (int i = 0; i < 4; i++){
        results.max_abs[i] = max_abs_diff[i];
        results.max_rel[i] = max_rel_diff[i];
    }

    if (std::getline(f1, line1) || std::getline(f2, line2)) {
        if (verbose) {
            report << "Row count mismatch detected\n\n";
        }
        anyMismatch = true;
    }
    is_equal = !anyMismatch;

    if (verbose) {
        report << "** Maximum Differences per Variable **\n";
        // Header
        report << std::left  << std::setw(12) << "Variable"
          << std::right << std::setw(18) << "MaxAbsDiff"
          << std::setw(18) << "MaxRelDiff" << "\n";
        report << std::string(48, '-') << "\n";
        report << std::setprecision(6) << std::scientific;
        for (int i=0; i<4; i++)
        {
            report << std::left << std::setw(12) << labels[i]
                << std::right << std::setw(18) << max_abs_diff[i] 
                << std::setw(18) << max_rel_diff[i] << "\n";
        }
        if (!is_equal)  report << RED << "==== ERROR: TEST NOT PASSED ====\n\n" << RESET;
        else report << GREEN << "==== TEST PASSED ====\n\n" << RESET;

        std::cout << report.str();
    }

    return results;
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

    // Flag for verbose output
    bool verbose = false;
	
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

        std::vector<CSVMaxDiff> file_summaries;
        int failed = 0;
		for (int i = 0; i < N_x; ++i)
		{
			for (int j = 0; j < N_y; ++j)
			{		
				rt_problem_ptr->write_surface_point_profiles_csv(test_out_path, i, j);

				if (mpi_rank == 0)
				{
					std::string test_file      = test_out_path + "_" + std::to_string(i) + "_" + std::to_string(j) + ".csv";
					std::string reference_file = (cfg.reference_sol_directory).string() + "profiles_" + std::to_string(i) + "_" + std::to_string(j) + ".csv";        
                    std::cout << "TEST: " << test_file << "\n";
                    std::cout << "REF : " << reference_file << "\n";
                    bool equal;
                    CSVMaxDiff report = compare_csv_report(test_file, reference_file, equal, verbose);
                    file_summaries.push_back(report);
                    if (!equal) failed++;
				}
			}
		}	
        
        if (failed > 0) {
            std::cout << RED << "\n==== ERROR: " << failed << "/" << N_x * N_y << " TESTS NOT PASSED ====\n" << RESET;
        } else {
            std::cout << GREEN << "\n==== TEST PASSED ====\n" << RESET;
        }


        if (mpi_rank == 0) {
            // Compute maximum over all files
            CSVMaxDiff global_max{};
            for (const auto& d : file_summaries) {
                for (int k=0; k<4; ++k) {
                    global_max.max_abs[k] = std::max(global_max.max_abs[k], d.max_abs[k]);
                    global_max.max_rel[k] = std::max(global_max.max_rel[k], d.max_rel[k]);
                }
            }
    
            // Print summary table
            std::cout << "\n** Maximum Differences Across All Files **\n";
            const char* labels[4] = {"I", "Q/I[%]", "U/I[%]", "V/I[%]"};
            // Header
            std::cout << std::left  << std::setw(12) << "Variable"
                    << std::right << std::setw(18) << "MaxAbsDiff"
                    << std::setw(18) << "MaxRelDiff" << "\n";
            std::cout << std::string(48, '-') << "\n";
            // Data rows
            std::cout << std::setprecision(6) << std::scientific;
            for (int i = 0; i < 4; ++i) {
                std::cout << std::left  << std::setw(12) << labels[i]
                        << std::right << std::setw(18) << global_max.max_abs[i]
                        << std::setw(18) << global_max.max_rel[i]
                        << "\n";
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
