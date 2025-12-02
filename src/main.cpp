#include "RT_solver.hpp"
#include "RT_utility.hpp"
#include "Test_rii_include.hpp"
#include "tools.h"

#include <chrono>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <string>
#include "tools.h"

// WARNING: if you want to use PORTA input for 3D setup, you need to set
// USE_PORTA_INPUT = 1 otherwise, it will use FAL-C input for 1D setup
#define USE_PORTA_INPUT 1
//////////////////////////////////////////////////////////////////////////

// WARNING: if you want to use command line options, you need to set
// USE_CMD_LINE_OPTIONS = 1 otherwise, it will use the default and hard-coded
// values
#define USE_CMD_LINE_OPTIONS 1
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
int
main(int argc, char *argv[])
{
	std::stringstream	  ss_a, ss_b;
	std::filesystem::path output_info_file_name;

	MPI_CHECK(MPI_Init(&argc, &argv));
	print_PETSc_mem();

	const double main_start_time = MPI_Wtime();

	int mpi_size;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	{ // check if the user wants to print the help message

		if (getOptionFlag(argc, argv, "--help"))
		{
			if (mpi_rank == 0) print_help();

			return 0;
		}
	}

	PetscInitialize(&argc, &argv, (char *)0, NULL);
	Kokkos::initialize(argc, argv);

#if ACC_SOLAR_3D == _ON_

#pragma message("ACC version ENABLED")

	const int devices_per_node = 4; // for MareNostrum 5 ACC & Daint Alps.
	int		  mpi_error		   = RII_epsilon_contrib::RII_contrib_MPI_Init(devices_per_node, MPI_COMM_WORLD);

	if (RII_epsilon_contrib::RII_contrib_MPI_Is_Device_Handler())
	{
		RII_epsilon_contrib::RII_contrib_MPI_Init_Memory_Pool();
	}

	for (int i = 0; i < LIMIT_OUT_DEVICE_MEMORY_USAGE; i++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank == i && RII_epsilon_contrib::RII_contrib_MPI_Is_Device_Handler())
		{
			std::cout << RII_epsilon_contrib::RII_contrib_MPI_Device_Info() << std::endl;
		}
	}

	const int is_device_handler = int(RII_epsilon_contrib::RII_contrib_MPI_Is_Device_Handler());

	int devices_cnt = 0;
	MPI_Reduce(&is_device_handler, &devices_cnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (mpi_rank == 0)
	{
		std::cout << "Using ACC version with " << devices_cnt << " devices." << std::endl;
	}

	set_RII_contrib_block_size(RII_CONTRIB_BLOCK_SIZE);

#else

	set_RII_contrib_block_size(1);

#endif // ACC_SOLAR_3D

	{ // start scope for RT_problem and RT_solver
		//////////////////////////////////////////////////////////////////////////
		// list of emissivity models
		// NONE:          undefined
		// CRD_limit:     CRD limit (Default set to CRD_GL)
		// CRD_limit_VHP: CRD limit with VHP (Very hight precision) approximation
		//
		// PRD:        default (PRD_FAST) partial redistribution. Default grid
		// (PRD_FAST). PRD_NORMAL: force partial redistribution with original grid.
		// PRD_FAST:   force partial redistribution with fast grid.
		//
		// PRD_AA:      partial redistribution Angle averaged method.
		// PRD_AA_MAPV: same as PRD_AA but it store the map of the values for a fast
		// computation
		//              (DANGER it uses a lot of memory)
		//
		// ZERO: continuum
		emissivity_model emissivity_model_var = emissivity_model::PRD;
		// emissivity_model emissivity_model_var = emissivity_model::PRD_NORMAL; //
		// Original grid over u'.

		const bool output					   = true;
		const bool output_overwrite_prevention = true; // if true the application stops (with an MPI_Abort) if
													   // the output directory already exists

		const bool use_B			= true;
		const bool use_ZERO_epsilon = false;

#if USE_CMD_LINE_OPTIONS == 1
		const bool use_CRD		 = getOptionFlag(argc, argv, "--CRD");
		const bool use_continuum = getOptionFlag(argc, argv, "--continuum");
		const bool use_prec		 = (not use_CRD);
#else
		const bool use_continuum = false;
		const bool use_CRD		 = true;
		const bool use_prec		 = (not use_CRD);
#endif

		if (use_CRD) emissivity_model_var = emissivity_model::CRD_limit;
		if (use_continuum) emissivity_model_var = emissivity_model::ZERO;

		// Set here the main input and output directories //////////////////////////
#if USE_CMD_LINE_OPTIONS == 1
		const std::filesystem::path main_input_dir	= getOptionArgument(argc, argv, "--input_dir");
		const std::filesystem::path main_output_dir = getOptionArgument(argc, argv, "--output_dir");
#else
		//  const std::filesystem::path main_input_dir  = "../input/PORTA";
		// const std::filesystem::path main_input_dir	= "/scratch/snx3000/pietrob/Comparison-TRIP-PORTA/a63/";
		// const std::filesystem::path main_output_dir = "/scratch/snx3000/pietrob/Comparison-TRIP-PORTA/a63/output";

		const std::filesystem::path main_input_dir	= "../input/PORTA/";
		const std::filesystem::path main_output_dir = "/gpfs/projects/ehpc238/output";
#endif
		////////////////////////////////////////////////////////////////////////////

#if USE_PORTA_INPUT == 1 // PORTA setup for 3D

// Set here the problem input file //////////////////////////
#if USE_CMD_LINE_OPTIONS == 1
		std::string input_pmd_string = getOptionArgument(argc, argv, "--problem_pmd_file");

		std::string input_cul_string;
		std::string input_qel_string;
		std::string input_llp_string;
		std::string input_back_string;

		const std::string input_config_string = getOptionArgument(argc, argv, "--problem_input_config");

		if (not input_config_string.empty())
		{
			const auto config_map =
				get_input_PORTA(main_input_dir / std::filesystem::path(input_config_string), mpi_rank);

			input_pmd_string = config_map.at("pmd");

			input_cul_string  = config_map.at("cul");
			input_qel_string  = config_map.at("qel");
			input_llp_string  = config_map.at("llp");
			input_back_string = config_map.at("back");
		}

#else
		const std::string input_pmd_string = std::string("AR_385_Cut_64x64_mirrorxy-CRD_I_V0-V0_fix_conv_KQ_MC.pmd");
		const std::string input_llp_string = std::string("AR_385_Cut_64x64_mirrorxy-CRD_I_V0-V0_fix_conv_KQ_MC.llp");

		const std::string input_cul_string	= std::string("AR_385_Cut_64x64_mirrorxy-CRD_I_V0_fix.cul");
		const std::string input_qel_string	= std::string("AR_385_Cut_64x64_mirrorxy-CRD_I_V0_fix.qel");
		const std::string input_back_string = std::string("AR_385_Cut_64x64_mirrorxy-CRD_I_V0_fix.back");
#endif
		/////////////////////////////////////////////////////////////

		// const auto input_pmd_file = std::filesystem::path(input_pmd_string);

		auto frequencies_input_path = main_input_dir / std::filesystem::path("frequency/96F");

		auto PORTA_input_pmd = main_input_dir / std::filesystem::path(input_pmd_string);

		// lambda to build the RT_problem object
		auto create_rt_problem = [&]()
		{
			if (input_cul_string.empty() or input_qel_string.empty() or input_llp_string.empty())
			{
				// VEECHIO  // solo PMD input at least one of the cul, qel, llp is
				// missing
				if (mpi_rank == 0)
				{
					std::cout << "WARNING: using ONLY PMD input file" << std::endl;
					std::cout << "Input PMD file: " << PORTA_input_pmd << std::endl;
				}

				return std::make_shared<RT_problem>(PORTA_input_pmd.string().c_str(), frequencies_input_path.string(),
													emissivity_model_var, use_B);
			}
			else
			{
				if (mpi_rank == 0)
					std::cout << "INPUT reading: using PMD + CUL + QEL + LLP + BACK input files" << std::endl;

				// NUOVO // PMD + CUL + QEL + LLP input
				// create cul and qel input path
				auto input_cul_path	 = main_input_dir / std::filesystem::path(input_cul_string);
				auto input_qel_path	 = main_input_dir / std::filesystem::path(input_qel_string);
				auto input_llp_path	 = main_input_dir / std::filesystem::path(input_llp_string);
				auto input_back_path = main_input_dir / std::filesystem::path(input_back_string);

				if (mpi_rank == 0)
				{
					const auto file_name_provided = [](const std::string &s)
					{ return (s.empty()) ? std::string("not provided") : s; };

					ss_a << "PMD input file:   " << file_name_provided(PORTA_input_pmd.string()) << std::endl;
					ss_a << "LLP input file:   " << file_name_provided(input_llp_path.string()) << std::endl;
					ss_a << "CUL input file:   " << file_name_provided(input_cul_path.string()) << std::endl;
					ss_a << "QEL input file:   " << file_name_provided(input_qel_path.string()) << std::endl;
					ss_a << "BACK input file:  " << file_name_provided(input_back_path.string()) << std::endl;
					std::cout << ss_a.str();
				}

				return std::make_shared<RT_problem>(PORTA_input_pmd.string().c_str(), input_cul_path.string().c_str(),
													input_qel_path.string().c_str(), input_llp_path.string().c_str(),
													input_back_path.string().c_str(), frequencies_input_path.string(),
													emissivity_model_var, use_B);
			}
		}; // end lambda create_rt_problem

		auto rt_problem_ptr = create_rt_problem(); /// CALL LAMBDA to create RT_problem object

		const int N_theta = rt_problem_ptr->N_theta_;
		const int N_chi	  = rt_problem_ptr->N_chi_;

		if (rt_problem_ptr->mpi_rank_ == 0)
		{
			// print the command line arguments
			ss_a << "MPI size = " << mpi_size << std::endl;
			ss_a << "Date and time: " << getCurrentDateTime() << std::endl;
			ss_a << std::endl << "Command line arguments: " << std::endl;
			ss_a << "argc: " << argc << std::endl;
			ss_a << "argv: ";
			for (int i = 0; i < argc; ++i)
			{
				ss_a << argv[i] << " ";
			}
			ss_a << std::endl << std::endl;

			ss_a << "PORTA 3D input file: " << PORTA_input_pmd << std::endl;
			ss_a << "Frequencies input path: " << frequencies_input_path << std::endl << std::endl;

			ss_a << "N_theta =            " << N_theta << std::endl;
			ss_a << "N_chi =              " << N_chi << std::endl << std::endl;

			const auto petsc_index_size = sizeof(PetscInt);
			ss_a << "PetscInt size: " << petsc_index_size << " bytes; " << (petsc_index_size * 8) << " bits." << std::endl
				 << std::endl;

			std::cout << ss_a.str();
		}

#else
		// FAL-C input for 1D input setup

		// inputs
		auto problem_input_file = std::filesystem::path("FAL-C/B20_V0_12T_8C_99F_1Pi4_9Pi8");
		auto problem_input_FAL	= main_input_dir / problem_input_file;
		// const std::string FAL_input_path = "../input/FAL-C/1_B0_V0_12T_8C_64F";
		// const std::string FAL_input_path =
		// "/users/pietrob/solar_3d/input/FAL-C/96F"; const std::string
		// FAL_input_path = "/users/pietrob/solar_3d/input/FAL-C/64F";

		const int N_theta	= 8;
		const int N_chi		= 16;
		auto rt_problem_ptr = std::make_shared<RT_problem>(problem_input_FAL.string(), N_theta, N_chi, use_CRD, use_B);
#endif

		RT_solver rt_solver(rt_problem_ptr, "BESSER", use_prec);
		// RT_solver rt_solver(rt_problem_ptr, "DELO_linear", use_prec);

		//////////////////////////////////////////////////////////////////////////
		// Prepare output directory
		// If the output directory does not exist, create it
		// It it exists, abort !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		std::string					output_file;
		const std::filesystem::path output_path =
			main_output_dir
			/ std::filesystem::path(input_pmd_string + "." + emissivity_model_to_string_long(emissivity_model_var));

		if (output) // Set up output
		{
			//  if (rt_problem_ptr->mpi_rank_ == 0)
			//    std::cout << "Output path: " << output_path << std::endl;

			if (rt_problem_ptr->mpi_rank_ == 0)
			{
				if (not std::filesystem::exists(output_path))
				{
					std::filesystem::create_directories(output_path);
				}
				else if (output_overwrite_prevention)
				{
					std::cerr << "File: " << __FILE__ << " Line: " << __LINE__
							  << " MPI rank: " << rt_problem_ptr->mpi_rank_ << " Directory: " << output_path
							  << " already exists" << std::endl;
					std::cerr << "Use different output directory" << std::endl;
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
			}

			// const std::string output_path =
			// "output/surface_profiles_64x64x133/B0V0/"; // TODO change
			output_file = (use_CRD) ? (output_path / "profiles_CRD").string() : (output_path / "profiles_PRD").string();

			if (rt_problem_ptr->mpi_rank_ == 0)
			{
				ss_b << "Output directory: " << output_path << std::endl;
				std::cout << ss_b.str();
				std::cout.flush();

				output_info_file_name = output_path / "info.txt"; // <-- Remove 'const auto', assign to outer variable
				std::ofstream output_file_info(output_info_file_name);

				output_file_info << ss_a.str();
				output_file_info << ss_b.str();
				output_file_info.close();
			}
		} // end if (output)

		const double main_setup_time = MPI_Wtime();
		if (rt_problem_ptr->mpi_rank_ == 0)
		{
			std::stringstream ss_mem;
			ss_mem << std::fixed << std::setprecision(2);
			ss_mem << "Setup time: " << (main_setup_time - main_start_time) << " seconds." << std::endl;
			ss_mem << "Total time: " << (main_setup_time - main_start_time) << " seconds." << std::endl;
			ss_mem << std::endl;

			std::cout << ss_mem.str();
			if (output) // Setp time output
			{
				std::ofstream output_file_info(output_info_file_name, std::ios_base::app);
				output_file_info << ss_mem.str();
				output_file_info.close();
			}
		}

		///////////////////////////////////////////////////
		// solve //////////////////////////////////////////
		rt_solver.solve();
		// rt_solver.apply_formal_solver();

		const double main_solve_end_time = MPI_Wtime();

		// lambda to compute arbitrary beam
		const auto compute_arbitrary_beam = [&, Nx = rt_problem_ptr->N_x_, Ny = rt_problem_ptr->N_y_](
												const Real mu, const Real chi, const std::string output_file)
		{
			char mu_charv[40];
			char chi_charv[40];

			std::snprintf(&mu_charv[0], 40, "%.4f", mu);
			std::snprintf(&chi_charv[0], 40, "%.4f", chi);

			std::string mu_str(mu_charv);
			std::string chi_str(chi_charv);

			// remove the dot from the output file name
			mu_str.erase(std::remove(mu_str.begin(), mu_str.end(), '.'), mu_str.end());
			chi_str.erase(std::remove(chi_str.begin(), chi_str.end(), '.'), chi_str.end());

			const std::string output_file_Omega_mu = output_file + "_mu" + mu_str + "_chi" + chi_str;
			// const std::string output_file_frequencies  = output_file + "_frequencies";
			// const std::string output_file_angular_grid = output_file + "_angular_grid";

			const Real theta = std::acos(mu);
			rt_solver.apply_formal_solver_Omega(theta, chi);

			for (int i = 0; i < Nx; ++i)
			{
				for (int j = 0; j < Ny; ++j)
				{
					rt_problem_ptr->write_surface_point_profiles_Omega(output_file_Omega_mu, i, j);
					rt_problem_ptr->write_surface_point_profiles_Omega_csv(output_file_Omega_mu, i, j, 14);
				}
			}
		}; // end lambda compute_arbitrary_beam

		// write output
		if (output) // Surface profiles for all surface points
		{
			const int N_x = rt_problem_ptr->N_x_;
			const int N_y = rt_problem_ptr->N_y_;

			for (int i = 0; i < N_x; ++i)
			{
				for (int j = 0; j < N_y; ++j)
				{
					// const std::string output_file_Omega_mu	   = output_file + "_mu" + mu_str + "_chi" + chi_str;
					const std::string output_file_frequencies  = (output_path / "frequencies_grid_Hz").string();
					const std::string output_file_angular_grid = (output_path / "angular_grid").string();

					rt_problem_ptr->write_surface_point_profiles(output_file, i, j);

					// TESTING purpose: write the CSV files too
					rt_problem_ptr->write_surface_point_profiles_csv(output_file, i, j);

					if (rt_problem_ptr->mpi_rank_ == 0 and i == 0 and j == 0)
					{
						rt_problem_ptr->write_angular_grid_csv(output_file_angular_grid, i, j);
						rt_problem_ptr->write_frequencies_grid_csv(output_file_frequencies, i, j);
					}
				}
			}
		} // end if (output)

		// rt_problem_ptr->write_surface_point_profiles(output_file, 0, 0);
		// rt_problem_ptr->write_surface_point_profiles(output_file, 10, 10);
		// rt_problem_ptr->write_surface_point_profiles(output_file, 3, 3);
		// rt_problem_ptr->write_surface_point_profiles(output_file, 1, 2);
		// rt_problem_ptr->write_surface_point_profiles(output_file, 1, 10);

		// old code: copied below .....

		// compute arbitrary beams
		////////////////////////////////////////////////////////////////////////////////////////////////
		/// Set arbitrary beam directions
		// free some memory
		rt_problem_ptr->free_fields_memory();
		rt_solver.free_fields_memory();

		// std::vector<Real>		mus_vec = {0.1, 0.3, 0.7, 1.0}; //// ATTENTION: arbitrary beam directions
		// const std::vector<Real> chi_vec = {0.0, 1.963495408493621e-01};

		std::vector<Real>		mus_vec = {0.1, 0.3, 0.7, 1.0}; //// ATTENTION: arbitrary beam directions
		const std::vector<Real> chi_vec = {0.0};

		if (rt_problem_ptr->mpi_rank_ == 0 and mus_vec.size() == 0)
		{
			std::cout << "WARNING: no arbitrary beams" << std::endl;
		}
		else if (rt_problem_ptr->mpi_rank_ == 0)
		{
// #define MU_EXTRA
#ifdef MU_EXTRA // extra beams
#pragma message "ATTENTION: using hardcoded extra arbitrary beams"
			{
				// const std::vector<Real> mus_extra = {0.183434642495650,
				// 0.960289856497537, 0.99}; const std::vector<Real> mus_extra =
				// {0.98, 0.99, 0.997};
				const std::vector<Real> mus_extra = {0.12, 0.3, 0.77}; // Test with different mu

				for (auto mux : mus_extra) mus_vec.push_back(mux);
				std::sort(mus_vec.begin(), mus_vec.end());
			}
#endif

			std::stringstream ss_ad;

			ss_ad << "\n" << std::string(50, '=') << std::endl;
			ss_ad << "Arbitrary Beam Directions" << std::endl;
			ss_ad << std::string(50, '=') << std::endl;
			ss_ad << std::setw(8) << "Beam #" << std::setw(15) << "mu" << std::setw(15) << "chi" << std::setw(12)
				  << "theta (rad)" << std::endl;
			ss_ad << std::string(50, '-') << std::endl;

			int beam_count = 0;
			for (Real chi : chi_vec)
			{
				for (Real mu : mus_vec)
				{
					Real theta = std::acos(mu);
					ss_ad << std::setw(8) << beam_count++ << std::setw(15) << std::fixed << std::setprecision(6) << mu
						  << std::setw(15) << std::fixed << std::setprecision(6) << chi << std::setw(12) << std::fixed
						  << std::setprecision(6) << theta << std::endl;
				}
			}
			ss_ad << std::string(50, '-') << std::endl;
			ss_ad << "Total number of beams: " << beam_count << std::endl;
			ss_ad << std::string(50, '=') << std::endl << std::endl;

			std::cout << ss_ad.str();

			if (!output_info_file_name.empty()) // <-- Add this check
			{
				std::ofstream output_file_info_ad(output_info_file_name, std::ios::app);
				output_file_info_ad << ss_ad.str();
				output_file_info_ad.close();
			} // end if !output_info_file_name.empty()

		} // end if mpi_rank_ == 0

		std::cout.flush();

		const double tick		   = MPI_Wtime();
		int			 beams_counter = 0;

		for (Real chi : chi_vec)
		{
			for (Real mu : mus_vec)
			{
				if (rt_problem_ptr->mpi_rank_ == 0)
				{
					std::cout << getCurrentDateTime() << " - Computing arbitrary beam: mu = " << mu << ", chi = " << chi
							  << std::endl;
				}

				compute_arbitrary_beam(mu, chi, output_file);
				beams_counter++;
			}
		}

		const double tock		   = MPI_Wtime();
		const double main_end_time = MPI_Wtime();

		if (rt_problem_ptr->mpi_rank_ == 0)
		{
			std::cout << "Arbitrary beams time (s) = " << tock - tick << std::endl;
			std::cout << "Number of beams          = " << beams_counter << std::endl;
			std::cout << "Time per beam (s)        = " << (tock - tick) / double(beams_counter) << std::endl;
		}

		if (output) // Final output and memory report
		{
#if ACC_SOLAR_3D == _ON_
			if (RII_epsilon_contrib::RII_contrib_MPI_Is_Device_Handler() and //
				mpi_rank < LIMIT_OUT_DEVICE_MEMORY_USAGE)					 //
			{
				std::string pool_info = RII_epsilon_contrib::RII_contrib_MPI_Generate_DPool_Report(true, true);

				const auto file_name_pool = output_path /													  //
											(boost::format("device_pool_info_rank_%d.txt") % mpi_rank).str(); //

				std::ofstream myfile;
				myfile.open(file_name_pool.string());
				if (myfile.is_open())
				{
					myfile << pool_info;
					myfile.close();
				}
				else
				{
					std::cerr << "Unable to open file: " << file_name_pool << std::endl;
				}
			}
#endif // ACC_SOLAR_3D

			// print memory usage
			const double byte_to_GB = 1.0 / (1000 * 1024 * 1024);

			unsigned long long vm_usage;
			unsigned long long resident_set;

			rii::process_mem_usage(vm_usage, resident_set);

			print_PETSc_mem();

			if (rt_problem_ptr->mpi_rank_ == 0)
			{
				std::stringstream ss_mem;
				ss_mem << "Total memory usage (vm_usage):               " << byte_to_GB * vm_usage << " GB" << std::endl;
				ss_mem << "Total memory usage (resident_set):           " << byte_to_GB * resident_set << " GB"
					   << std::endl;

				print_PETSc_mem();
				// ss_mem << mem_petsc << std::endl << std::endl;

#if ACC_SOLAR_3D == _ON_
				ss_mem << "Total number of devices (accelerators) used: " << devices_cnt << std::endl;
#endif // ACC_SOLAR_3D
				ss_mem << "MPI size:                                    " << mpi_size << std::endl;

				ss_mem << std::fixed << std::setprecision(2) << std::endl;
				ss_mem << std::endl;
				ss_mem << "╔════════════════════════════════════════════════╗" << std::endl;
				ss_mem << "║           EXECUTION TIME SUMMARY               ║" << std::endl;
				ss_mem << "╠════════════════════════════════════════════════╣" << std::endl;
				ss_mem << "║ Setup time:            " << std::setw(10) << (main_setup_time - main_start_time)
					   << " seconds      ║" << std::endl;
				ss_mem << "║ Solve time:            " << std::setw(10) << (main_solve_end_time - main_setup_time)
					   << " seconds      ║" << std::endl;
				ss_mem << "║ Post processing time:  " << std::setw(10) << (main_end_time - main_solve_end_time)
					   << " seconds      ║" << std::endl;
				ss_mem << "╠════════════════════════════════════════════════╣" << std::endl;
				ss_mem << "║ Total execution time:  " << std::setw(10) << (main_end_time - main_start_time)
					   << " seconds      ║" << std::endl;
				ss_mem << "╚════════════════════════════════════════════════╝" << std::endl;
				ss_mem << std::endl;

				std::cout << ss_mem.str();

				if (!output_info_file_name.empty()) // <-- Add this check
				{
					std::ofstream output_file_info(output_info_file_name, std::ios::app);
					output_file_info << ss_mem.str();
					output_file_info.close();
				} // end if !output_info_file_name.empty()
			} // end if mpi_rank_ == 0

		} // end if (output)
	} // end scope for RT_problem and RT_solver

#if ACC_SOLAR_3D == _ON_
	if (RII_epsilon_contrib::RII_contrib_MPI_Is_Device_Handler())
		RII_epsilon_contrib::RII_contrib_MPI_Finalize_Memory_Pool();

	RII_epsilon_contrib::RII_contrib_MPI_Finalize();
#endif // ACC_SOLAR_3D

	Kokkos::finalize();
	PetscFinalize(); // CHKERRQ(ierr);
	MPI_CHECK(MPI_Finalize());

	return 0;
} // end main
//////////////////////////////////////////////////////////////////////////